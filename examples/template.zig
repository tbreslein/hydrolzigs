// Copyright (c) 2023-2024
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

const std = @import("std");
const hydrol = @import("hydrol");

const conf = hydrol.Config{
    .mesh = .{
        .type = .cartesian,
        .n = 10,
        .xi_in = 1.0,
        .xi_out = 2.0,
    },
    .physics = .{
        .type = .euler1d_isothermal,
    },
    .numflux = .{
        .limiter_mode = .vanleer,
    },
};

pub fn main() !void {
    const mesh = hydrol.Mesh(conf){};
    var u = hydrol.Physics(conf){};
    var rhs = hydrol.RHS(conf){};
    hydrol.run(conf, &u, &rhs, mesh);
}
