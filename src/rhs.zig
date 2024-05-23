// Copyright (c) 2023-2024
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the RHS struct, responsible for solving the right hand side of the equations.

const Config = @import("config.zig").Config;
const getNumEq = @import("physics.zig").getNumEq;
const Mesh = @import("mesh.zig").Mesh;
const vsplat = @import("utils.zig").vsplat;
const NumFlux = @import("numflux.zig").NumFlux;
const Physics = @import("physics.zig").Physics;

/// Handles solving the right hand side of the set of differential equations.
///
/// Takes a Config as a comptime value which stores the configuration for a simulation.
pub fn RHS(comptime c: Config) type {
    const num_eq = getNumEq(c);
    return struct {
        full_rhs: [num_eq][c.mesh.n]f64 = [_]@Vector(c.mesh.n, f64){@splat(0.0)} ** num_eq,
        numflux: NumFlux(c) = .{},

        pub fn updateRHS(self: *RHS(c), u: *Physics(c), mesh: Mesh(c)) void {
            self.*.numflux.calcDFluxDXi(u);
            inline for (0..num_eq) |j| {
                inline for (mesh.ixi_in..mesh.ixi_out + 1) |i| {
                    self.*.full_rhs[j][i] = self.*.numflux.flux_num[j][i];
                }
            }
        }
    };
}
