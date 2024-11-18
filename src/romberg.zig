fn Romberg(KMAXX: usize) type {
    if (KMAXX < 2) @compileError("KMAXX should be at least 2.");

    const slice_al: comptime_int = @alignOf([]f64);
    const child_al: comptime_int = @alignOf(f64);
    const slice_sz: comptime_int = @sizeOf(usize) * 2;
    const child_sz: comptime_int = @sizeOf(f64);
    const coeff_sz: comptime_int = ((KMAXX * KMAXX + KMAXX) >> 1) * child_sz + KMAXX * slice_sz;

    return struct {
        coeff: [][]f64,
        table: []f64,

        const Self = @This();

        const Method = enum { once, iter };

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

        fn integrate(self: *const Self, comptime method: Method, f: *const fn (x: f64) f64, a: f64, b: f64) f64 {
            switch (method) {
                .once => {
                    for (0..KMAXX) |n| self.table[n] = trapezoidal_sum(f, a, b, (n + 1) << 1);

                    const offset: f64 = self.table[KMAXX - 1];
                    for (0..KMAXX) |n| self.table[n] -= offset;

                    for (1..KMAXX) |k| {
                        for (0..KMAXX - k, 1..) |i, ip1| {
                            self.table[i] += self.coeff[k][i] * (self.table[i] - self.table[ip1]);
                        }
                    }

                    return offset + self.table[0];
                },
                .iter => {
                    self.table[0] = trapezoidal_sum(f, a, b, (0 + 1) << 1);
                    for (1..KMAXX) |n| {
                        self.table[n] = trapezoidal_sum(f, a, b, (n + 1) << 1);
                        var k: usize = n - 1;
                        while (true) : (k -= 1) {
                            self.table[k] += self.coeff[n - k][k] * (self.table[k] - self.table[k + 1]);
                            if (k == 0) break;
                        }
                    }

                    return self.table[0];
                },
            }
        }
    };
}

test "Romberg" {
    const page = testing.allocator;

    const self = try Romberg(8).init(page);
    defer self.deinit(page);

    debug.print("{d}\n", .{self.integrate(.once, func, 0.0, 2.5)});
    debug.print("{d}\n", .{self.integrate(.iter, func, 0.0, 2.5)});
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
