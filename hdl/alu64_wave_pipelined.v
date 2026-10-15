`timescale 1ps/100fs

// 64-Bit Wave-Pipelined Kogge-Stone RSFQ Arithmetic Logic Unit (ALU)
// Reference: Inamdar et al., IEEE Trans. Appl. Supercond., Vol. 31, No. 5, August 2021
// Total Junctions: 90,256 (64-Bit), 6,882 (8-Bit Reference)
// Architecture: Wave-pipelined Kogge-Stone Carry Tree (INIT, STR & REG, SUM)

// -----------------------------------------------------------------------------
// Sub-Block 1: Kogge-Stone Prefix Black Cell & Gray Cell
// -----------------------------------------------------------------------------
module ksa_black_cell #(
    parameter real T_PD = 6.5
) (
    input  wire g_high,
    input  wire p_high,
    input  wire g_low,
    input  wire p_low,
    output wire g_out,
    output wire p_out
);
    // G_out = G_high | (P_high & G_low)
    // P_out = P_high & P_low
    assign #(T_PD) g_out = g_high | (p_high & g_low);
    assign #(T_PD) p_out = p_high & p_low;
endmodule

module ksa_gray_cell #(
    parameter real T_PD = 6.0
) (
    input  wire g_high,
    input  wire p_high,
    input  wire g_low,
    output wire g_out
);
    // G_out = G_high | (P_high & G_low)
    assign #(T_PD) g_out = g_high | (p_high & g_low);
endmodule

// -----------------------------------------------------------------------------
// Sub-Block 2: INIT Sub-Block (Operand Pre-conditioning & Logic/Bitwise Generation)
// -----------------------------------------------------------------------------
module alu_init_subblock #(
    parameter WIDTH = 64,
    parameter real T_INIT = 8.5
) (
    input  wire             clk,
    input  wire [2:0]       opcode,   // 000: ADD, 001: SUB, 010: AND, 011: OR, 100: XOR, 101: NOT
    input  wire [WIDTH-1:0] a,
    input  wire [WIDTH-1:0] b,
    output reg  [WIDTH-1:0] g,        // Bitwise Generate: A & B
    output reg  [WIDTH-1:0] p,        // Bitwise Propagate: A | B
    output reg  [WIDTH-1:0] p_sum,    // Bitwise Partial Sum: A ^ B
    output reg  [WIDTH-1:0] logic_out // Direct logic output for logic mode
);
    integer i;
    reg [WIDTH-1:0] b_eff;

    always @(posedge clk) begin
        #(T_INIT);
        for (i = 0; i < WIDTH; i = i + 1) begin
            // Operand conditioning (Invert B for subtraction)
            b_eff[i] = (opcode == 3'b001) ? ~b[i] : b[i];

            // Bitwise Generate, Propagate, and Partial Sum
            g[i]     <= a[i] & b_eff[i];
            p[i]     <= a[i] | b_eff[i];
            p_sum[i] <= a[i] ^ b_eff[i];

            // Logic operations evaluated in INIT stage
            case (opcode)
                3'b010: logic_out[i] <= a[i] & b[i]; // AND
                3'b011: logic_out[i] <= a[i] | b[i]; // OR
                3'b100: logic_out[i] <= a[i] ^ b[i]; // XOR
                3'b101: logic_out[i] <= ~a[i];       // NOT
                default: logic_out[i] <= a[i] ^ b_eff[i];
            endcase
        end
    end
endmodule

// -----------------------------------------------------------------------------
// Sub-Block 3: 64-Bit Kogge-Stone Parallel-Prefix Network (STR & REG Stages)
// -----------------------------------------------------------------------------
module alu_kogge_stone_prefix_tree_64 #(
    parameter real T_STAGE = 7.0
) (
    input  wire        clk,
    input  wire [63:0] g_in,
    input  wire [63:0] p_in,
    input  wire        cin,
    output reg  [63:0] carry_out
);
    // 6 Prefix Levels: log2(64) = 6
    reg [63:0] g_l1, p_l1;
    reg [63:0] g_l2, p_l2;
    reg [63:0] g_l3, p_l3;
    reg [63:0] g_l4, p_l4;
    reg [63:0] g_l5, p_l5;
    reg [63:0] g_l6, p_l6;

    integer i;

    // Level 1: Step = 1
    always @(*) begin
        for (i = 0; i < 64; i = i + 1) begin
            if (i >= 1) begin
                g_l1[i] = g_in[i] | (p_in[i] & g_in[i-1]);
                p_l1[i] = p_in[i] & p_in[i-1];
            end else begin
                g_l1[0] = g_in[0] | (p_in[0] & cin);
                p_l1[0] = p_in[0];
            end
        end
    end

    // Level 2: Step = 2
    always @(*) begin
        for (i = 0; i < 64; i = i + 1) begin
            if (i >= 2) begin
                g_l2[i] = g_l1[i] | (p_l1[i] & g_l1[i-2]);
                p_l2[i] = p_l1[i] & p_l1[i-2];
            end else begin
                g_l2[i] = g_l1[i];
                p_l2[i] = p_l1[i];
            end
        end
    end

    // Level 3: Step = 4
    always @(*) begin
        for (i = 0; i < 64; i = i + 1) begin
            if (i >= 4) begin
                g_l3[i] = g_l2[i] | (p_l2[i] & g_l2[i-4]);
                p_l3[i] = p_l2[i] & p_l2[i-4];
            end else begin
                g_l3[i] = g_l2[i];
                p_l3[i] = p_l2[i];
            end
        end
    end

    // Level 4: Step = 8
    always @(*) begin
        for (i = 0; i < 64; i = i + 1) begin
            if (i >= 8) begin
                g_l4[i] = g_l3[i] | (p_l3[i] & g_l3[i-8]);
                p_l4[i] = p_l3[i] & p_l3[i-8];
            end else begin
                g_l4[i] = g_l3[i];
                p_l4[i] = p_l3[i];
            end
        end
    end

    // Level 5: Step = 16
    always @(*) begin
        for (i = 0; i < 64; i = i + 1) begin
            if (i >= 16) begin
                g_l5[i] = g_l4[i] | (p_l4[i] & g_l4[i-16]);
                p_l5[i] = p_l4[i] & p_l4[i-16];
            end else begin
                g_l5[i] = g_l4[i];
                p_l5[i] = p_l4[i];
            end
        end
    end

    // Level 6: Step = 32 (Final level for 64 bits)
    always @(*) begin
        for (i = 0; i < 64; i = i + 1) begin
            if (i >= 32) begin
                g_l6[i] = g_l5[i] | (p_l5[i] & g_l5[i-32]);
                p_l6[i] = p_l5[i] & p_l5[i-32];
            end else begin
                g_l6[i] = g_l5[i];
                p_l6[i] = p_l5[i];
            end
        end
    end

    always @(posedge clk) begin
        #(6 * T_STAGE);
        carry_out <= g_l6;
    end
endmodule

// -----------------------------------------------------------------------------
// Sub-Block 4: Top-Level 64-Bit Wave-Pipelined RSFQ ALU
// -----------------------------------------------------------------------------
module alu64_top #(
    parameter real T_CLK_PERIOD = 40.0 // Default 40 ps (25 GHz Nominal)
) (
    input  wire        clk,
    input  wire [2:0]  opcode,
    input  wire [63:0] a,
    input  wire [63:0] b,
    input  wire        cin,
    output reg  [63:0] result,
    output reg         cout
);
    wire [63:0] g_init, p_init, p_sum, logic_out;
    wire [63:0] carries;

    // Instantiate INIT Sub-Block
    alu_init_subblock #(.WIDTH(64)) u_init (
        .clk(clk),
        .opcode(opcode),
        .a(a),
        .b(b),
        .g(g_init),
        .p(p_init),
        .p_sum(p_sum),
        .logic_out(logic_out)
    );

    // Instantiate 64-Bit Kogge-Stone Tree
    alu_kogge_stone_prefix_tree_64 u_ksa (
        .clk(clk),
        .g_in(g_init),
        .p_in(p_init),
        .cin(cin),
        .carry_out(carries)
    );

    // SUM Sub-Block: Final output multiplexing between Arithmetic & Logic
    integer i;
    always @(posedge clk) begin
        if (opcode == 3'b000 || opcode == 3'b001) begin
            // Arithmetic mode: ADD / SUB
            result[0] <= p_sum[0] ^ cin;
            for (i = 1; i < 64; i = i + 1) begin
                result[i] <= p_sum[i] ^ carries[i-1];
            end
            cout <= carries[63];
        end else begin
            // Logic mode: AND, OR, XOR, NOT
            result <= logic_out;
            cout <= 1'b0;
        end
    end
endmodule
