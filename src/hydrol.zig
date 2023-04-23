// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Entry point for the hydrolzig library.
//! Reexports config, as well as the init and run functions.

const std = @import("std");
pub const Config = @import("config.zig").Config;
pub const Mesh = @import("mesh.zig").Mesh;
pub const Physics = @import("physics.zig").Physics;

/// Run a simulation
pub fn run(comptime c: Config, u: *Physics(c), mesh: Mesh(c)) void {
    std.debug.print("mesh = {}\n", .{mesh});
    std.debug.print("u = {}\n", .{u});
}
