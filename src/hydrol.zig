const std = @import("std");
pub const config = @import("config.zig");
const m_mesh = @import("mesh.zig");

pub fn init(comptime conf: config.Config) m_mesh.Mesh(conf.mesh.n) {
    comptime conf.validate();
    return m_mesh.init_mesh(conf);
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
