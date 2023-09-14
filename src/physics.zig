// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the Physics struct that carries the state of physics variables

const Config = @import("config.zig").Config;
const Mesh = @import("mesh.zig").Mesh;
const vsplat = @import("mesh.zig").vsplat;
const std = @import("std");
const panic = std.debug.panic;
const math = std.math;

/// Returns the number of a equations for a simulation configuration stored in the comptime parameter c.
pub fn getNumEq(comptime c: Config) u32 {
    return switch (c.physics.type) {
        .euler1d_adiabatic => 3,
        .euler1d_isothermal => 2,
    };
}

/// Handles the variables regarding the physics in the simulation.
///
/// Takes a Config as a comptime value which stores the configuration for a simulation.
pub fn Physics(comptime c: Config) type {
    const num_eq = getNumEq(c);
    const j_rho = switch (c.physics.type) {
        .euler1d_adiabatic => 0,
        .euler1d_isothermal => 0,
    };
    const j_xi = switch (c.physics.type) {
        .euler1d_adiabatic => 1,
        .euler1d_isothermal => 1,
    };
    const j_pressure = switch (c.physics.type) {
        .euler1d_adiabatic => 2,
        .euler1d_isothermal => 2,
    };
    const j_eigenmin = 0;
    const j_eigenmax = num_eq - 1;

    const is_isothermal = switch (c.physics.type) {
        .euler1d_adiabatic => false,
        .euler1d_isothermal => true,
    };

    return struct {
        /// Variables at the centre of each cell in the mesh
        cent: Variables = .{},

        /// Variables at the west faces of each cell in the mesh
        west: Variables = .{},

        /// Variables at the east faces of each cell in the mesh
        east: Variables = .{},

        /// Equation index for the minimal eigen values in Variables.eigen_vals
        comptime j_eigenmin: u32 = j_eigenmin,

        /// Equation index for the maximal eigen values in Variables.eigen_vals
        comptime j_eigenmax: u32 = j_eigenmax,

        pub const Variables = struct {
            /// Primitive variables
            prim: [num_eq]@Vector(c.mesh.n, f64) = [_]@Vector(c.mesh.n, f64){@splat(0.0)} ** num_eq,

            /// Conservative variables
            cons: [num_eq]@Vector(c.mesh.n, f64) = [_]@Vector(c.mesh.n, f64){@splat(0.0)} ** num_eq,

            /// Speed of sound
            csound: @Vector(c.mesh.n, f64) = @splat(0.0),

            /// Eigen values
            eigen_vals: [num_eq]@Vector(c.mesh.n, f64) = [_]@Vector(c.mesh.n, f64){@splat(0.0)} ** num_eq,

            /// Physical flux
            flux: [num_eq]@Vector(c.mesh.n, f64) = [_]@Vector(c.mesh.n, f64){@splat(0.0)} ** num_eq,

            /// Assigns an instance of Variables to self
            pub fn assignFrom(self: *Variables, other: Variables) void {
                self.*.prim = other.prim;
                self.*.cons = other.cons;
                self.*.csound = other.csound;
            }

            /// Update conservative variables, assuming the primitive ones are up-to-date
            pub fn updateCons(self: *Variables) void {
                switch (c.physics.type) {
                    .euler1d_isothermal => {
                        self.*.cons[j_rho] = self.*.prim[j_rho];
                        self.*.cons[j_xi] = self.*.prim[j_xi] * self.*.prim[j_rho];
                    },
                    .euler1d_adiabatic => {
                        self.*.cons[j_rho] = self.*.prim[j_rho];
                        self.*.cons[j_xi] = self.*.prim[j_xi] * self.*.prim[j_rho];
                        self.*.cons[j_pressure] = self.*.prim[j_pressure] / vsplat(c.mesh.n, c.physics.adiabatic_index - 1.0) + (vsplat(c.mesh.n, 0.5) * self.*.prim[j_rho] * self.*.prim[j_xi] * self.*.prim[j_xi]);
                    },
                }
            }

            /// Update primitive variables, assuming the conservative ones are up-to-date
            pub fn updatePrim(self: *Variables) void {
                switch (c.physics.type) {
                    .euler1d_isothermal => {
                        self.*.prim[j_rho] = self.*.cons[j_rho];
                        self.*.prim[j_xi] = self.*.cons[j_xi] / self.*.cons[j_rho];
                    },
                    .euler1d_adiabatic => {
                        self.*.prim[j_rho] = self.*.cons[j_rho];
                        self.*.prim[j_xi] = self.*.cons[j_xi] / self.*.cons[j_rho];
                        self.*.prim[j_pressure] = vsplat(c.mesh.n, c.physics.adiabatic_index - 1.0) * (self.*.cons[j_pressure] - vsplat(c.mesh.n, 0.5) / self.*.cons[j_rho] * self.*.cons[j_xi] * self.*.cons[j_xi]);
                    },
                }
            }

            /// Update speed of sound, assuming the primitive variables are up-to-date
            pub fn updateCsound(self: *Variables) void {
                if (!is_isothermal) {
                    self.*.csound = @sqrt(self.*.gamma_vec * self.*.prim[j_pressure] / self.*.prim[j_rho]);
                }
            }

            /// Update eigen values, assuming the primitive variables are up-to-date
            pub fn updateEigenVals(self: *Variables) void {
                self.*.eigen_vals[j_eigenmin] = self.*.prim[j_xi] - self.*.csound;
                self.*.eigen_vals[j_eigenmax] = self.*.prim[j_xi] + self.*.csound;
                for (self.*.eigen_vals[1..j_eigenmax]) |*ev| {
                    ev.* = self.*.prim[j_xi];
                }
            }

            /// Update physical flux, assuming the primitive variables are up-to-date
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

        /// Initialises self.west and self.east, assuming self.cent is already initialised.
        ///
        /// This is especially useful in isothermal physics, because updateCsound is a no-op in that case,
        /// so even if self.cent.csound is set, those fields in self.{west,east} need to be set through this
        /// function because it would never be set otherwise.
        pub fn initWestEast(self: *Physics(c)) void {
            self.*.west.assignFrom(self.*.cent);
            self.*.east.assignFrom(self.*.cent);
        }

        /// Calculates the time step width according to the CFL criterium.
        ///
        /// Panics in case the time step width is not finite.
        pub fn calcDtCFL(self: Physics(c), c_cfl: f64, mesh: Mesh) f64 {
            const dt = c_cfl / @reduce(.Max, @fabs(self.cent.eigen_vals[j_eigenmax] * mesh.cell_width_inv));
            if (!math.isFinite(dt)) {
                panic("CFL time step turned non-finite!");
            }
            return dt;
        }
    };
}

fn conversionTest(comptime physconf: Config.PhysicsConfig) !void {
    const c = comptime Config{
        .mesh = .{
            .type = .cartesian,
            .n = 10,
            .xi_in = 1.0,
            .xi_out = 2.0,
        },
        .physics = physconf,
        .numflux = .{
            .limiter_mode = .vanleer,
        },
    };
    const m = comptime getNumEq(c);
    var expected = [_]@Vector(c.mesh.n, f64){vsplat(c.mesh.n, 111.1)} ** m;
    var u = Physics(c){};
    u.cent.prim = expected;
    u.cent.updateCons();
    u.cent.updatePrim();
    for (0..m) |j| {
        for (0..c.mesh.n) |i| {
            try std.testing.expectApproxEqRel(expected[j][i], u.cent.prim[j][i], 1.0E-12);
        }
    }
}

test "Convert back and forth // euler1d isothermal" {
    try conversionTest(.{ .type = .euler1d_isothermal });
}

test "Convert back and forth // euler1d adiabatic" {
    try conversionTest(.{ .type = .euler1d_adiabatic, .adiabatic_index = 5.0 / 3.0 });
}
