fn Romberg(KMAXX: usize) type {
    if (KMAXX < 2) @compileError("KMAXX should be at least 2.");
    if (KMAXX & 1 == 1) @compileError("KMAXX should be even.");

    const slice_al: comptime_int = @alignOf([]f64);
    const child_al: comptime_int = @alignOf(f64);
    const slice_sz: comptime_int = @sizeOf(usize) * 2;
    const child_sz: comptime_int = @sizeOf(f64);
    const coeff_sz: comptime_int = (KMAXX + 1) * (KMAXX >> 1) * child_sz + KMAXX * slice_sz;

    return struct {
        coeff: [][]f64,
        table: []f64,

        const Self = @This();

        fn init(allocator: mem.Allocator) !*Self {
            const self: *Self = try allocator.create(Self);
            errdefer allocator.destroy(self);

            self.table = try allocator.alloc(f64, KMAXX);
            errdefer allocator.free(self.table);

            self.coeff = outer: {
                const buff: []u8 = try allocator.alloc(u8, coeff_sz);
                errdefer allocator.free(buff);

                const temp: [][]f64 = inner: {
                    const ptr: [*]align(slice_al) []f64 = @ptrCast(@alignCast(buff.ptr));
                    break :inner ptr[0..KMAXX];
                };

                var padding: usize = comptime KMAXX * slice_sz;
                var chunk_sz: usize = comptime KMAXX * child_sz;

                for (0..KMAXX) |k| {
                    temp[k] = inner: {
                        const ptr: [*]align(child_al) f64 = @ptrCast(@alignCast(buff.ptr + padding));
                        break :inner ptr[0 .. KMAXX - k];
                    };
                    padding += chunk_sz;
                    chunk_sz -= child_sz;
                }

                break :outer temp;
            };

            for (self.coeff[0], 0..) |*p, i| p.* = @as(f64, @floatFromInt((i + 1) << 1));
            for (self.coeff[1..], 1..) |row, k| {
                for (row, 0.., k..) |*p, i, ipk| {
                    p.* = 1.0 / (sqr(self.coeff[0][i] / self.coeff[0][ipk]) - 1.0);
                }
            }
            return self;
        }

        fn deinit(self: *const Self, allocator: mem.Allocator) void {
            {
                const ptr: [*]u8 = @ptrCast(@alignCast(self.coeff.ptr));
                const len: usize = coeff_sz;
                allocator.free(ptr[0..len]);
            }

            allocator.free(self.table);
            allocator.destroy(self);
        }

        fn integrate(self: *const Self, f: *const fn (x: f64) f64, a: f64, b: f64) f64 {
            for (self.table, 0..) |*p, n| p.* = trapezoidal_sum(f, a, b, (n + 1) << 1);

            const offset: f64 = self.table[KMAXX - 1];
            for (self.table) |*p| p.* -= offset;

            for (1..KMAXX) |k| {
                for (0..KMAXX - k, 1..) |i, ip1| {
                    self.table[i] += self.coeff[k][i] * (self.table[i] - self.table[ip1]);
                }
            }

            return offset + self.table[0];
        }
    };
}

test "Romberg" {
    const page = testing.allocator;

    const self = try Romberg(4).init(page);
    defer self.deinit(page);

    for (self.coeff) |row| debug.print("{any}\n", .{row});

    const result: f64 = self.integrate(func, 0.0, 2.5);
    debug.print("{d}\n", .{result});
}

fn trapezoidal_sum(f: *const fn (x: f64) f64, a: f64, b: f64, n: usize) f64 {
    const step: f64 = (b - a) / @as(f64, @floatFromInt(n));
    const half: f64 = 0.5 * step;

    var x_lo: f64 = a;
    var x_hi: f64 = undefined;
    var f_lo: f64 = undefined;
    var f_hi: f64 = undefined;

    var sum: f64 = 0.0;

    for (0..n) |_| {
        x_hi = x_lo + step;
        f_lo = f(x_lo);
        f_hi = f(x_hi);
        sum += half * (f_lo + f_hi);
        x_lo = x_hi;
    }

    return sum;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

fn sqr(x: f64) f64 {
    return x * x;
}

fn func(x: f64) f64 {
    return (x + 1) * (x - 1) * (x - 2);
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

const std = @import("std");
const mem = std.mem;
const math = std.math;
const debug = std.debug;
const testing = std.testing;
