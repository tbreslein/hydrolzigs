// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the RHS struct, responsible for solving the right hand side of the equations.

const Config = @import("config.zig").Config;
const getNumEq = @import("physics.zig").getNumEq;
const Mesh = @import("mesh.zig").Mesh;
const NumFlux = @import("numflux.zig").NumFlux;
const Physics = @import("physics.zig").Physics;

/// Handles solving the right hand side of the set of differential equations.
///
/// Takes a Config as a comptime value which stores the configuration for a simulation.
pub fn RHS(comptime c: Config) type {
    const num_eq = getNumEq(c);
    return struct {
        full_rhs: [num_eq]@Vector(c.mesh.n, f64) = [_]@Vector(c.mesh.n, f64){@splat(c.mesh.n, @as(f64, 0.0))} ** num_eq,
        numflux: NumFlux(c) = .{},
    };
}
