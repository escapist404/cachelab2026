/*
请注意，你的代码不能出现任何 int/short/char/float/double/auto 等局部变量/函数传参，我们仅允许使用 reg 定义的寄存器变量。
其中 reg 等价于一个 int。

你不能自己申请额外的内存，即不能使用 new/malloc，作为补偿我们传入了一段 buffer，大小为 BUFFER_SIZE = 64，你可以视情况使用。

我们的数组按照 A, B, C, buffer 的顺序在内存上连续紧密排列，且 &A = 0x30000000（这是模拟的设定，不是 A 的真实地址）

如果你需要以更自由的方式访问内存，你可以以相对 A 的方式访问，比如 A[100]，用 *(0x30000000) 是无法访问到的。

如果你有定义常量的需求（更严谨的说法是，你想定义的是汇编层面的立即数，不应该占用寄存器），请参考下面这种方式使用宏定义来完成。
*/

#include "cachelab.h"

#define m case0_m
#define n case0_n
#define p case0_p

// 我们用这个 2*2*2 的矩阵乘法来演示寄存器是怎么被分配的
void gemm_case0(ptr_reg A, ptr_reg B, ptr_reg C, ptr_reg buffer) {  // allocate 0 1 2 3
    for (reg i = 0; i < m; ++i) {                                   // allocate 4
        for (reg j = 0; j < p; ++j) {                               // allocate 5
            reg tmpc = 0;                                           // allocate 6
            for (reg k = 0; k < n; ++k) {                           // allocate 7
                reg tmpa = A[i * n + k];                            // allocate 8
                reg tmpb = B[k * p + j];                            // allocate 9
                tmpc += tmpa * tmpb;
            }  // free 9 8
            // free 7
            C[i * p + j] = tmpc;
        }  // free 6
        // free 5
    }
    // free 4
}  // free 3 2 1 0

#undef m
#undef n
#undef p

#define m case1_m
#define n case1_n
#define p case1_p

#define CALC_COL_4(c_arr, cols, col)     \
    b = B[k * p + j + (col)];            \
    (c_arr)[0 * (cols) + (col)] += a[0] * b; \
    (c_arr)[1 * (cols) + (col)] += a[1] * b; \
    (c_arr)[2 * (cols) + (col)] += a[2] * b; \
    (c_arr)[3 * (cols) + (col)] += a[3] * b

void gemm_case1(ptr_reg A, ptr_reg B, ptr_reg C, ptr_reg buffer) {
    for (reg i = 0; i < m; i += 4) {
        for (reg j_step = 0; j_step < p; j_step += 4) {
            reg j = j_step;
            if ((i / 4) % 2 != 0) {
                j = p - 4 - j_step;
            }
            reg c[16]{};

            if (((i / 4) * 4 + (j_step / 4)) % 2 == 0) {
                for (reg k = 0; k < n; ++k) {
                    reg a[4];
                    reg b[4];
                    for (reg ii = 0; ii < 4; ++ii) {
                        a[ii] = A[(i + ii) * n + k];
                    }
                    for (reg ii = 0; ii < 4; ++ii) {
                        b[ii] = B[k * p + j + ii];
                    }
                    for (reg ii = 0; ii < 4; ++ii) {
                        for (reg jj = 0; jj < 4; ++jj) {
                            c[ii * 4 + jj] += a[ii] * b[jj];
                        }
                    }
                }
            } else {
                for (reg k = n - 1; k >= 0; --k) {
                    reg a[4];
                    reg b[4];
                    for (reg ii = 0; ii < 4; ++ii) {
                        a[ii] = A[(i + ii) * n + k];
                    }
                    for (reg ii = 0; ii < 4; ++ii) {
                        b[ii] = B[k * p + j + ii];
                    }
                    for (reg ii = 0; ii < 4; ++ii) {
                        for (reg jj = 0; jj < 4; ++jj) {
                            c[ii * 4 + jj] += a[ii] * b[jj];
                        }
                    }
                }
            }

            for (reg ii = 0; ii < 4; ++ii) {
                for (reg jj = 0; jj < 4; ++jj) {
                    C[(i + ii) * p + j + jj] = c[ii * 4 + jj];
                }
            }
        }
    }
}

#undef m
#undef n
#undef p

#define m case2_m
#define n case2_n
#define p case2_p

void gemm_case2(ptr_reg A, ptr_reg B, ptr_reg C, ptr_reg buffer) {
    for (reg i = 0; i < m; i += 4) {
        for (reg j = 0; j < p; j += 4) {
            reg c[16]{};

            if ((j / 4) % 2 == 0) {
                for (reg k = 0; k < n; ++k) {
                    reg a[4];
                    reg b[4];
                    for (reg ii = 0; ii < 4; ++ii) {
                        a[ii] = A[(i + ii) * n + k];
                    }
                    for (reg ii = 0; ii < 4; ++ii) {
                        b[ii] = B[k * p + j + ii];
                    }
                    for (reg ii = 0; ii < 4; ++ii) {
                        for (reg jj = 0; jj < 4; ++jj) {
                            c[ii * 4 + jj] += a[ii] * b[jj];
                        }
                    }
                }
            } else {
                for (reg k = n - 1; k >= 0; --k) {
                    reg a[4];
                    reg b[4];
                    for (reg ii = 0; ii < 4; ++ii) {
                        a[ii] = A[(i + ii) * n + k];
                    }
                    for (reg ii = 0; ii < 4; ++ii) {
                        b[ii] = B[k * p + j + ii];
                    }
                    for (reg ii = 0; ii < 4; ++ii) {
                        for (reg jj = 0; jj < 4; ++jj) {
                            c[ii * 4 + jj] += a[ii] * b[jj];
                        }
                    }
                }
            }

            for (reg ii = 0; ii < 4; ++ii) {
                for (reg jj = 0; jj < 4; ++jj) {
                    C[(i + ii) * p + j + jj] = c[ii * 4 + jj];
                }
            }
        }
    }
}

#undef m
#undef n
#undef p

#define m case3_m
#define n case3_n
#define p case3_p

void gemm_case3(ptr_reg A, ptr_reg B, ptr_reg C, ptr_reg buffer) {
    for (reg j = 0; j < p - 3; j += 6) {
        for (reg i = 0; i < m - 3; i += 4) {
            reg c[24]{};

            if ((i / 4) % 2 == 0) {
                for (reg k = 0; k < n; ++k) {
                    reg a[4];
                    a[0] = A[i * n + k];
                    a[1] = A[(i + 1) * n + k];
                    a[2] = A[(i + 2) * n + k];
                    a[3] = A[(i + 3) * n + k];

                    reg b;
                    CALC_COL_4(c, 6, 0);
                    CALC_COL_4(c, 6, 1);
                    CALC_COL_4(c, 6, 2);
                    CALC_COL_4(c, 6, 3);
                    CALC_COL_4(c, 6, 4);
                    CALC_COL_4(c, 6, 5);
                }
            } else {
                for (reg k = n - 1; k >= 0; --k) {
                    reg a[4];
                    a[0] = A[i * n + k];
                    a[1] = A[(i + 1) * n + k];
                    a[2] = A[(i + 2) * n + k];
                    a[3] = A[(i + 3) * n + k];

                    reg b;
                    CALC_COL_4(c, 6, 0);
                    CALC_COL_4(c, 6, 1);
                    CALC_COL_4(c, 6, 2);
                    CALC_COL_4(c, 6, 3);
                    CALC_COL_4(c, 6, 4);
                    CALC_COL_4(c, 6, 5);
                }
            }

            for (reg ii = 0; ii < 4; ++ii) {
                for (reg jj = 0; jj < 6; ++jj) {
                    C[(i + ii) * p + j + jj] = c[ii * 6 + jj];
                }
            }
        }

        {
            reg i = 24;
            reg c[18]{};

            for (reg k = 0; k < n; ++k) {
                reg a[3];
                reg b[6];
                for (reg ii = 0; ii < 3; ++ii) {
                    a[ii] = A[(i + ii) * n + k];
                }
                for (reg jj = 0; jj < 6; ++jj) {
                    b[jj] = B[k * p + j + jj];
                }
                for (reg ii = 0; ii < 3; ++ii) {
                    for (reg jj = 0; jj < 6; ++jj) {
                        c[ii * 6 + jj] += a[ii] * b[jj];
                    }
                }
            }

            for (reg ii = 0; ii < 3; ++ii) {
                for (reg jj = 0; jj < 6; ++jj) {
                    C[(i + ii) * p + j + jj] = c[ii * 6 + jj];
                }
            }
        }
    }

    {
        reg j = 24;
        for (reg i = 0; i < m - 3; i += 4) {
            reg c[12]{};

            if ((i / 4) % 2 == 0) {
                for (reg k = 0; k < n; ++k) {
                    reg a[4];
                    reg b[3];
                    for (reg ii = 0; ii < 4; ++ii) {
                        a[ii] = A[(i + ii) * n + k];
                    }
                    for (reg jj = 0; jj < 3; ++jj) {
                        b[jj] = B[k * p + j + jj];
                    }
                    for (reg ii = 0; ii < 4; ++ii) {
                        for (reg jj = 0; jj < 3; ++jj) {
                            c[ii * 3 + jj] += a[ii] * b[jj];
                        }
                    }
                }
            } else {
                for (reg k = n - 1; k >= 0; --k) {
                    reg a[4];
                    reg b[3];
                    for (reg ii = 0; ii < 4; ++ii) {
                        a[ii] = A[(i + ii) * n + k];
                    }
                    for (reg jj = 0; jj < 3; ++jj) {
                        b[jj] = B[k * p + j + jj];
                    }
                    for (reg ii = 0; ii < 4; ++ii) {
                        for (reg jj = 0; jj < 3; ++jj) {
                            c[ii * 3 + jj] += a[ii] * b[jj];
                        }
                    }
                }
            }

            for (reg ii = 0; ii < 4; ++ii) {
                for (reg jj = 0; jj < 3; ++jj) {
                    C[(i + ii) * p + j + jj] = c[ii * 3 + jj];
                }
            }
        }

        {
            reg i = 24;
            reg c[9]{};

            for (reg k = 0; k < n; ++k) {
                reg a[3];
                reg b[3];
                for (reg ii = 0; ii < 3; ++ii) {
                    a[ii] = A[(i + ii) * n + k];
                }
                for (reg jj = 0; jj < 3; ++jj) {
                    b[jj] = B[k * p + j + jj];
                }
                for (reg ii = 0; ii < 3; ++ii) {
                    for (reg jj = 0; jj < 3; ++jj) {
                        c[ii * 3 + jj] += a[ii] * b[jj];
                    }
                }
            }

            for (reg ii = 0; ii < 3; ++ii) {
                for (reg jj = 0; jj < 3; ++jj) {
                    C[(i + ii) * p + j + jj] = c[ii * 3 + jj];
                }
            }
        }
    }
}

#undef m
#undef n
#undef p 