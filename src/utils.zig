// Copyright (c) 2023-2024
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports utility functions

pub fn vsplat(comptime n: usize, comptime x: f64) @Vector(n, f64) {
    return @splat(x);
}
