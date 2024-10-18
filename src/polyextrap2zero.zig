fn polyExtrapZero(allocator: mem.Allocator, xs: []f64, ys: []f64, n: usize) !f64 {
    const offset: f64 = blk: {
        var min: f64 = math.inf(f64);
        var tmp: f64 = undefined;
        var ind: usize = undefined;
        for (0..n) |i| {
            tmp = @abs(0.0 - xs[i]);
            if (tmp == 0.0) return ys[i];
            if (tmp < min) {
                min = tmp;
                ind = i;
            }
        }
        break :blk ys[ind];
    };

    const table: []f64 = try allocator.alloc(f64, n);
    defer allocator.free(table);

    for (0..n) |i| table[i] = ys[i] - offset;

    for (1..n) |j| {
        for (0..n - j, 1..) |i, ip1| {
            table[i] += (table[i] - table[ip1]) / (xs[i + j] / xs[i] - 1.0);
        }
    }

    return offset + table[0];
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

fn func(x: f64) f64 {
    return (x + 1) * (x - 1) * (x - 2);
}

test "polyExtrapZero" {
    const page = testing.allocator;
    var xs: [4]f64 = .{ -1.5, 0.5, 1.5, 2.5 };
    var ys: [4]f64 = undefined;
    for (xs, &ys) |x, *y| y.* = func(x);
    try testing.expectApproxEqAbs(func(0.0), try polyExtrapZero(page, &xs, &ys, 4), 3e-16);
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

const std = @import("std");
const mem = std.mem;
const math = std.math;
const testing = std.testing;

const Array = @import("./array.zig").Array;
