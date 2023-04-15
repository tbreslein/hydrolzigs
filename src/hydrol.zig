//! Entry point for the hydrolzig library.
//! Reexports config, as well as the init and run functions.

const std = @import("std");
pub const config = @import("config.zig");
const m_mesh = @import("mesh.zig");

/// Initialises all objects needed for a simulation, configured by a hydrol.config.Config
pub fn init(comptime conf: config.Config) m_mesh.Mesh(conf.mesh.n) {
    comptime conf.validate();
    return m_mesh.initMesh(conf);
}

/// Run a simulation, given the objects you initialised with init()
pub fn run(comptime N: u32, _: m_mesh.Mesh(N)) void {}
