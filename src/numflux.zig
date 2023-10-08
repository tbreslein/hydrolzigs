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
    const n = c.mesh.n;
    const n_reduced = n - 1;
    const m = getNumEq(c);
    const mesh = Mesh(c){};
    const inv_dxi: f64 = 1.0 / mesh.dxi;
    const dist_west: @Vector(n, f64) = mesh.xi_west - mesh.xi_cent;
    const dist_east: @Vector(n, f64) = mesh.xi_east - mesh.xi_cent;
    const da_east = blk: {
        var vec = vsplat(n_reduced, 0.0);
        inline for (0.., mesh.ixi_in - 1..mesh.ixi_out + 1) |i, j| {
            vec[i] = mesh.d_area_xi_deta_dphi_east[j];
        }
        break :blk vec;
    };

    return struct {
        slopes: @Vector(n, f64) = vsplat(n, 0.0),
        a_plus: @Vector(n_reduced, f64) = vsplat(n_reduced, 0.0),
        a_minus: @Vector(n_reduced, f64) = vsplat(n_reduced, 0.0),
        b: @Vector(n_reduced, f64) = vsplat(n_reduced, 0.0),
        c: @Vector(n_reduced, f64) = vsplat(n_reduced, 0.0),
        flux_num: [m]@Vector(n_reduced, f64) = [_]@Vector(n_reduced, f64){vsplat(n_reduced, 0.0)} ** m,
        eigen_west: [m]@Vector(n_reduced, f64) = [_]@Vector(n_reduced, f64){vsplat(n_reduced, 0.0)} ** m,
        eigen_east: [m]@Vector(n_reduced, f64) = [_]@Vector(n_reduced, f64){vsplat(n_reduced, 0.0)} ** m,
        uflux_west: [m]@Vector(n_reduced, f64) = [_]@Vector(n_reduced, f64){vsplat(n_reduced, 0.0)} ** m,
        uflux_east: [m]@Vector(n_reduced, f64) = [_]@Vector(n_reduced, f64){vsplat(n_reduced, 0.0)} ** m,
        ucons_west: [m]@Vector(n_reduced, f64) = [_]@Vector(n_reduced, f64){vsplat(n_reduced, 0.0)} ** m,
        ucons_east: [m]@Vector(n_reduced, f64) = [_]@Vector(n_reduced, f64){vsplat(n_reduced, 0.0)} ** m,

        pub fn calcDFluxDXi(self: *NumFlux(c), u: *Physics(c)) void {
            self.*.reconstruct(u);
            u.*.west.updatePrim();
            u.*.east.updatePrim();
            u.*.west.updateCSound();
            u.*.east.updateCSound();
            u.*.west.updateEigenVals();
            u.*.east.updateEigenVals();
            u.*.west.updateFlux();
            u.*.east.updateFlux();

            inline for (0..m) |j| {
                inline for (0.., mesh.ixi_in - 1..mesh.ixi_out + 1, mesh.ixi_in..mesh.ixi_out + 2) |i_0, i_1, i_2| {
                    self.*.eigen_east[j][i_0] = u.*.east.eigen_vals[j][i_1];
                    self.*.eigen_west[j][i_0] = u.*.west.eigen_vals[j][i_2];
                    self.*.uflux_east[j][i_0] = u.*.east.flux[j][i_1];
                    self.*.uflux_west[j][i_0] = u.*.west.flux[j][i_2];
                    self.*.ucons_east[j][i_0] = u.*.east.cons[j][i_1];
                    self.*.ucons_west[j][i_0] = u.*.west.cons[j][i_2];
                }
            }

            // TODO: store the eigen vals, a+, a-, b, c, fluxnum, uflux_east, uflux_west, ucons_east, ucons_west,
            // and d_area_xi_deta_dphi_east in fitting Vectors
            // I could just have vectors for these values that are all just n_comp + 1 in size, and load the values
            // I need directly into those
            // AND benchmark it against a regular for loop version

            self.*.a_plus = @max(vsplat(n_reduced, 0.0), @max(self.*.eigen_west[m - 1], self.eigen_east[m - 1]));
            self.*.a_minus = @min(vsplat(n_reduced, 0.0), @min(self.*.eigen_west[0], self.eigen_east[0]));
            self.*.b = da_east / (self.*.a_plus - self.*.a_minus);
            self.*.c = self.*.a_plus * self.*.a_minus;
        }

        fn reconstruct(self: *NumFlux(c), u: *Physics(c)) void {
            inline for (0..m) |j| {
                inline for (1..n - 1) |i| {
                    const x = u.*.cent.cons[j][i] - u.*.cent.cons[j][i - 1];
                    const y = u.*.cent.cons[j][i - 1] - u.*.cent.cons[j][i];

                    self.*.slopes[i] = inv_dxi * switch (c.numflux.limiter_mode) {
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
                            break :blk (x * abs_x + y * abs_y) / (abs_x + abs_y + math.floatMin(f64));
                        },
                    };
                }
                u.*.west.cons[j] = u.*.cent.cons[j] + self.*.slopes * dist_west;
                u.*.east.cons[j] = u.*.cent.cons[j] + self.*.slopes * dist_east;
            }
        }
    };
}
