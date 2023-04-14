const std = @import("std");
pub const config = @import("config.zig");

pub fn init(comptime conf: config.Config) u32 {
    comptime conf.validate();
    return conf.mesh.n_gc;
}

pub fn add(x: i32, y: i32) i32 {
    return x + y;
}

test "simple test" {
    var list = std.ArrayList(i32).init(std.testing.allocator);
    defer list.deinit(); // try commenting this out and see if zig detects the memory leak!
    try list.append(42);
    try std.testing.expectEqual(@as(i32, 42), list.pop());
}
