// Copyright (c) 2023-2024
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Entry point for the hydrolzig library.
//! Reexports config, as well as the init and run functions.

const std = @import("std");
pub const Config = @import("config.zig").Config;
pub const Mesh = @import("mesh.zig").Mesh;
pub const Physics = @import("physics.zig").Physics;
pub const RHS = @import("rhs.zig").RHS;

/// Run a simulation
pub fn run(comptime c: Config, u: *Physics(c), rhs: *RHS(c), mesh: Mesh(c)) void {
    std.debug.print("\n", .{});
    std.debug.print("mesh.type = {}\n\n", .{mesh.type});
    rhs.*.updateRHS(u, mesh);
    std.debug.print("flux_num = {any}\n\n", .{rhs.numflux.flux_num});
    std.debug.print("full_rhs = {any}\n\n", .{rhs.full_rhs});
}

test {
    @import("std").testing.refAllDecls(@This());
}
