// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the NumFlux struct, responsible for calculating the numerical flux

const Config = @import("config.zig").Config;
const getNumEq = @import("physics.zig").getNumEq;
const Mesh = @import("mesh.zig").Mesh;
const vsplat = @import("utils.zig").vsplat;
const Physics = @import("physics.zig").Physics;
const std = @import("std");
const math = std.math;

/// Handles calculating the numerical flux according to a generalised Kurganov-Tadmor scheme.
///
/// Takes a Config as a comptime value which stores the configuration for a simulation.
pub fn NumFlux(comptime c: Config) type {
    const num_eq = getNumEq(c);
    const mesh_temp = Mesh(c){};
    return struct {
        a_plus: @Vector(c.mesh.n, f64) = vsplat(c.mesh.n, 0.0),
        a_minus: @Vector(c.mesh.n, f64) = vsplat(c.mesh.n, 0.0),
        b: @Vector(c.mesh.n, f64) = vsplat(c.mesh.n, 0.0),
        c: @Vector(c.mesh.n, f64) = vsplat(c.mesh.n, 0.0),
        comptime inv_dxi: f64 = 1.0 / mesh_temp.dxi,
        comptime dist_west: @Vector(c.mesh.n, f64) = mesh_temp.xi_west - mesh_temp.xi_cent,
        comptime dist_east: @Vector(c.mesh.n, f64) = mesh_temp.xi_east - mesh_temp.xi_cent,
        flux_num: [num_eq]@Vector(c.mesh.n, f64) = [_]@Vector(c.mesh.n, f64){vsplat(c.mesh.n, 0.0)} ** num_eq,

        pub fn calcDFluxDXi(self: *NumFlux, u: *Physics, _: Mesh) void {
            self.*.reconstruct(u);
        }

        fn reconstruct(self: NumFlux, u: *Physics(c)) void {
            inline for (0..num_eq) |j| {
                inline for (1..c.mesh.n - 1) |i| {
                    const x = u.*.cent.cons[j][i] - u.*.cent.cons[j][i - 1];
                    const y = u.*.cent.cons[j][i - 1] - u.*.cent.cons[j][i];

                    // this should be resolved at comptime, since c is comptime
                    const slope = self.inv_dxi * comptime switch (c.numflux.limiter_mode) {
                        .minmod => if (math.sign(x) * math.sign(y) > 0) {
                            math.sign(x) * math.min(@fabs(x), @fabs(y));
                        } else {
                            0.0;
                        },
                        .superbee => if (x * y > 0) {
                            math.min(math.min(@fabs(x), @fabs(y)), 0.5 * math.max(@fabs(x), @fabs(y)));
                        } else {
                            0.0;
                        },
                        .monocent => blk: {
                            const z = 0.5 * (u.*.cent.cons[j][i - 1] - u.*.cent.cons[j][i + 1]);
                            break :blk if (math.sign(x) * math.sign(y) > 0 and math.sign(y) * math.sign(z) > 0) {
                                math.sign(x) * math.min(@fabs(x * c.numflux.limiter_param), math.min(@fabs(y * c.numflux.limiter_param), @fabs(y)));
                            } else {
                                0.0;
                            };
                        },
                        .vanleer => blk: {
                            const abs_x = @fabs(x);
                            const abs_y = @fabs(y);
                            break :blk (x * abs_x + y * abs_y) / (abs_x + abs_y + math.f64_min);
                        },
                    };
                    u.*.west.cons[j][i] = u.*.cent.cons[j][i] + slope * self.dist_west;
                    u.*.east.cons[j][i] = u.*.cent.cons[j][i] + slope * self.dist_east;
                }
            }
        }
    };
}
