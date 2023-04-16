const Config = @import("config.zig").Config;
const Mesh = @import("mesh.zig").Mesh;
const std = @import("std");
const panic = std.debug.panic;
const math = std.math;

pub fn Physics(comptime c: Config) type {
    const num_eq = switch (c.physics.type) {
        .euler1d_adiabatic => 3,
        .euler1d_isothermal => 2,
    };
    const j_rho = switch (c.physics.type) {
        .euler1d_adiabatic => 1,
        .euler1d_isothermal => 1,
    };
    const j_xi = switch (c.physics.type) {
        .euler1d_adiabatic => 2,
        .euler1d_isothermal => 2,
    };
    const j_pressure = switch (c.physics.type) {
        .euler1d_adiabatic => 3,
        .euler1d_isothermal => 3,
    };
    const j_eigenmin = 0;
    const j_eigenmax = num_eq - 1;

    const is_isothermal = switch (c.physics.type) {
        .euler1d_adiabatic => false,
        .euler1d_isothermal => true,
    };

    return struct {
        cent: Variables = .{},
        west: Variables = .{},
        east: Variables = .{},

        comptime j_eigenmin: u32 = j_eigenmin,
        comptime j_eigenmax: u32 = j_eigenmax,

        pub const Variables = struct {
            prim: [num_eq]@Vector(c.mesh.n, f64) = [_]@Vector(c.mesh.n, f64){@splat(c.mesh.n, @as(f64, 0.0))} ** num_eq,
            cons: [num_eq]@Vector(c.mesh.n, f64) = [_]@Vector(c.mesh.n, f64){@splat(c.mesh.n, @as(f64, 0.0))} ** num_eq,
            csound: @Vector(c.mesh.n, f64) = @splat(c.mesh.n, @as(f64, 0.0)),
            eigen_vals: [num_eq]@Vector(c.mesh.n, f64) = [_]@Vector(c.mesh.n, f64){@splat(c.mesh.n, @as(f64, 0.0))} ** num_eq,
            flux: [num_eq]@Vector(c.mesh.n, f64) = [_]@Vector(c.mesh.n, f64){@splat(c.mesh.n, @as(f64, 0.0))} ** num_eq,

            pub fn assignFrom(self: *Variables, other: Variables) void {
                self.*.prim = other.prim;
                self.*.cons = other.cons;
                self.*.csound = other.csound;
            }

            pub fn updateCons(self: *Variables) void {
                switch (c.physics.type) {
                    .euler1d_isothermal => {
                        self.*.cons[j_rho] = self.*.prim[j_rho];
                        self.*.cons[j_xi] = self.*.prim[j_xi] * self.*.prim[j_rho];
                    },
                    .euler1d_adiabatic => {
                        self.*.cons[j_rho] = self.*.prim[j_rho];
                        self.*.cons[j_xi] = self.*.prim[j_xi] * self.*.prim[j_rho];
                        self.*.cons[j_pressure] = self.*.prim[j_pressure] / @splat(c.mesh.n, c.physics.adiabatic_index - 1.0) + (0.5 * self.*.prim[j_rho] * self.*.prim[j_xi] * self.*.prim[j_xi]);
                    },
                }
            }

            pub fn updatePrim(self: *Variables) void {
                switch (c.physics.type) {
                    .euler1d_isothermal => {
                        self.*.prim[j_rho] = self.*.cons[j_rho];
                        self.*.prim[j_xi] = self.*.cons[j_xi] / self.*.cons[j_rho];
                    },
                    .euler1d_adiabatic => {
                        self.*.prim[j_rho] = self.*.cons[j_rho];
                        self.*.prim[j_xi] = self.*.cons[j_xi] / self.*.cons[j_rho];
                        self.*.prim[j_pressure] = @splat(c.mesh.n, c.physics.adiabatic_index - 1.0) * (self.*.cons[j_pressure] - @splat(c.mesh.n, 0.5) / self.*.cons[j_rho] * self.*.cons[j_xi] * self.*.cons[j_xi]);
                    },
                }
            }

            pub fn updateCsound(self: *Variables) void {
                if (!is_isothermal) {
                    self.*.csound = @sqrt(@splat(c.mesh.n, c.physics.adiabatic_index) * self.*.prim[j_pressure] / self.*.prim[j_rho]);
                }
            }

            pub fn updateEigenVals(self: *Variables) void {
                self.*.eigen_vals[j_eigenmin] = self.*.prim[j_xi] - self.*.csound;
                self.*.eigen_vals[j_eigenmax] = self.*.prim[j_xi] + self.*.csound;
                for (self.*.eigen_vals[1..j_eigenmax]) |*ev| {
                    ev.* = self.*.prim[j_xi];
                }
            }

            pub fn updateFlux(self: *Variables) void {
                switch (c.physics.type) {
                    .euler1d_isothermal => {
                        self.*.flux[j_rho] = self.*.cons[j_xi];
                        self.*.flux[j_xi] = self.*.cons[j_rho] * (self.*.prim[j_xi] * self.*.prim[j_xi] + self.*.csound * self.*.csound);
                    },
                    .euler1d_adiabatic => {
                        self.*.flux[j_rho] = self.*.cons[j_xi];
                        self.*.flux[j_xi] = @mulAdd(@Vector(c.mesh.n, f64), self.*.prim[j_xi], self.*.cons[j_xi], self.*.prim[j_pressure]);
                        self.*.flux[j_pressure] = self.*.prim[j_xi] * (self.*.cons[j_pressure] + self.*.prim[j_pressure]);
                    },
                }
            }
        };

        pub fn initWestEast(self: *Physics(c)) void {
            self.*.west.assignFrom(self.*.cent);
            self.*.east.assignFrom(self.*.cent);
        }

        pub fn calcDtCFL(self: Physics(c), c_cfl: f64, mesh: Mesh) f64 {
            const dt = c_cfl / @reduce(.Max, @fabs(self.cent.eigen_vals[j_eigenmax] * mesh.cell_width_inv));
            if (!math.isFinite(dt)) {
                panic("CFL time step turned non-finite!");
            }
            return dt;
        }
    };
}
