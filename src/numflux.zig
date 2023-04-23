// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the NumFlux struct, responsible for calculating the numerical flux

const Config = @import("config.zig").Config;
const getNumEq = @import("physics.zig").getNumEq;
const Mesh = @import("mesh.zig").Mesh;
const Physics = @import("physics.zig").Physics;

/// Handles calculating the numerical flux according to a generalised Kurganov-Tadmor scheme.
///
/// Takes a Config as a comptime value which stores the configuration for a simulation.
pub fn NumFlux(comptime c: Config) type {
    const num_eq = getNumEq(c);
    const mesh_temp = Mesh(c){};
    return struct {
        a_plus: @Vector(c.mesh.n, f64) = @splat(c.mesh.n, @as(f64, 0.0)),
        a_minus: @Vector(c.mesh.n, f64) = @splat(c.mesh.n, @as(f64, 0.0)),
        b: @Vector(c.mesh.n, f64) = @splat(c.mesh.n, @as(f64, 0.0)),
        c: @Vector(c.mesh.n, f64) = @splat(c.mesh.n, @as(f64, 0.0)),
        inv_dxi: f64 = 1.0 / mesh_temp.dxi,
        dist_west: @Vector(c.mesh.n, f64) = mesh_temp.xi_west - mesh_temp.xi_cent,
        dist_east: @Vector(c.mesh.n, f64) = mesh_temp.xi_east - mesh_temp.xi_cent,
        flux_num: [num_eq]@Vector(c.mesh.n, f64) = [_]@Vector(c.mesh.n, f64){@splat(c.mesh.n, @as(f64, 0.0))} ** num_eq,

        pub fn calcDFluxDXi(self: *NumFlux, u: *Physics, _: Mesh) void {
            self.*.reconstruct(u);
        }

        fn reconstruct(_: NumFlux, _: *Physics(c)) void {
            for (0..num_eq) |_| {}
        }
    };
}
