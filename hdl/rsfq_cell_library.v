`timescale 1ps/100fs

// Synthesizable RSFQ Standard Cell Library for MIT-LL SFQ5ee Process
// Reference: Inamdar et al., IEEE Trans. Appl. Supercond., Vol. 31, No. 5, August 2021
// Supports parameterizable propagation delays for Timing Back-Annotation

// -----------------------------------------------------------------------------
// 1. D-Flip-Flop (dff_00a) with Internal Storage Loop
// -----------------------------------------------------------------------------
module rsfq_dff #(
    parameter real T_CLK_TO_Q = 8.5,
    parameter real T_SETUP    = 3.0,
    parameter real T_HOLD     = 2.0
) (
    input  wire clk,
    input  wire d,
    output reg  q
);
    reg state; // 0: Empty, 1: Flux Stored

    initial begin
        state = 0;
        q = 0;
    end

    always @(posedge d) begin
        state <= 1'b1;
    end

    always @(posedge clk) begin
        if (state == 1'b1) begin
            #(T_CLK_TO_Q) q <= 1'b1;
            #2.0 q <= 1'b0; // Single-fluxon pulse output (width ~2 ps)
            state <= 1'b0;  // Destructive readout resets loop
        end
    end
endmodule

// -----------------------------------------------------------------------------
// 2. D-Flip-Flop with Complementary Outputs (dffc_01a)
// -----------------------------------------------------------------------------
module rsfq_dffc #(
    parameter real T_CLK_TO_Q  = 8.5,
    parameter real T_CLK_TO_QI = 8.8
) (
    input  wire clk,
    input  wire d,
    output reg  q,
    output reg  qi
);
    reg state;

    initial begin
        state = 0;
        q = 0;
        qi = 0;
    end

    always @(posedge d) begin
        state <= 1'b1;
    end

    always @(posedge clk) begin
        if (state == 1'b1) begin
            #(T_CLK_TO_Q) q <= 1'b1;
            #2.0 q <= 1'b0;
            state <= 1'b0;
        end else begin
            #(T_CLK_TO_QI) qi <= 1'b1;
            #2.0 qi <= 1'b0;
        end
    end
endmodule

// -----------------------------------------------------------------------------
// 3. Inverter / NOT Cell (not_00c)
// -----------------------------------------------------------------------------
module rsfq_not #(
    parameter real T_CLK_TO_OUT = 7.5
) (
    input  wire clk,
    input  wire in,
    output reg  out
);
    reg state; // 1: Primed to emit on clk, 0: Inhibited by input pulse

    initial begin
        state = 1'b1;
        out = 0;
    end

    always @(posedge in) begin
        state <= 1'b0;
    end

    always @(posedge clk) begin
        if (state == 1'b1) begin
            #(T_CLK_TO_OUT) out <= 1'b1;
            #2.0 out <= 1'b0;
        end
        state <= 1'b1; // Re-arm for next cycle
    end
endmodule

// -----------------------------------------------------------------------------
// 4. Clocked AND Cell (and_00c)
// -----------------------------------------------------------------------------
module rsfq_and #(
    parameter real T_CLK_TO_OUT = 9.2
) (
    input  wire clk,
    input  wire a,
    input  wire b,
    output reg  out
);
    reg state_a, state_b;

    initial begin
        state_a = 0;
        state_b = 0;
        out = 0;
    end

    always @(posedge a) state_a <= 1'b1;
    always @(posedge b) state_b <= 1'b1;

    always @(posedge clk) begin
        if (state_a && state_b) begin
            #(T_CLK_TO_OUT) out <= 1'b1;
            #2.0 out <= 1'b0;
        end
        state_a <= 1'b0;
        state_b <= 1'b0;
    end
endmodule

// -----------------------------------------------------------------------------
// 5. Destructive / Data AND Cell (dand_00c)
// -----------------------------------------------------------------------------
module rsfq_dand #(
    parameter real T_PD = 6.5
) (
    input  wire a,
    input  wire b,
    output reg  out
);
    reg state_b;

    initial begin
        state_b = 0;
        out = 0;
    end

    always @(posedge b) state_b <= 1'b1;

    always @(posedge a) begin
        if (state_b) begin
            #(T_PD) out <= 1'b1;
            #2.0 out <= 1'b0;
            state_b <= 1'b0;
        end
    end
endmodule

// -----------------------------------------------------------------------------
// 6. Clocked XOR Cell (xor_00a)
// -----------------------------------------------------------------------------
module rsfq_xor #(
    parameter real T_CLK_TO_OUT = 8.8
) (
    input  wire clk,
    input  wire a,
    input  wire b,
    output reg  out
);
    reg state;

    initial begin
        state = 0;
        out = 0;
    end

    always @(posedge a) state <= ~state;
    always @(posedge b) state <= ~state;

    always @(posedge clk) begin
        if (state == 1'b1) begin
            #(T_CLK_TO_OUT) out <= 1'b1;
            #2.0 out <= 1'b0;
            state <= 1'b0;
        end
    end
endmodule

// -----------------------------------------------------------------------------
// 7. Demultiplexer / 1-to-2 Router (d2_02a)
// -----------------------------------------------------------------------------
module rsfq_d2 #(
    parameter real T_PD = 7.0
) (
    input  wire clk,
    input  wire d,
    output reg  q0,
    output reg  q1
);
    reg state;

    initial begin
        state = 0;
        q0 = 0;
        q1 = 0;
    end

    always @(posedge d) state <= 1'b1;

    always @(posedge clk) begin
        if (state == 1'b1) begin
            #(T_PD) q1 <= 1'b1;
            #2.0 q1 <= 1'b0;
            state <= 1'b0;
        end else begin
            #(T_PD) q0 <= 1'b1;
            #2.0 q0 <= 1'b0;
        end
    end
endmodule

// -----------------------------------------------------------------------------
// 8. Active Splitter (s22a_00a)
// -----------------------------------------------------------------------------
module rsfq_splitter #(
    parameter real T_PD = 3.5
) (
    input  wire in,
    output reg  out1,
    output reg  out2
);
    initial begin
        out1 = 0;
        out2 = 0;
    end

    always @(posedge in) begin
        #(T_PD) begin
            out1 <= 1'b1;
            out2 <= 1'b1;
        end
        #2.0 begin
            out1 <= 1'b0;
            out2 <= 1'b0;
        end
    end
endmodule

// -----------------------------------------------------------------------------
// 9. Josephson Transmission Lines (JTL1 & JTL2)
// -----------------------------------------------------------------------------
module rsfq_jtl1 #(
    parameter real T_PD = 2.2
) (
    input  wire in,
    output reg  out
);
    initial out = 0;
    always @(posedge in) begin
        #(T_PD) out <= 1'b1;
        #2.0 out <= 1'b0;
    end
endmodule

module rsfq_jtl2 #(
    parameter real T_PD = 4.4
) (
    input  wire in,
    output reg  out
);
    initial out = 0;
    always @(posedge in) begin
        #(T_PD) out <= 1'b1;
        #2.0 out <= 1'b0;
    end
endmodule

// -----------------------------------------------------------------------------
// 10. Half-Adder Cell (ha_03a)
// -----------------------------------------------------------------------------
module rsfq_ha #(
    parameter real T_SUM = 8.5,
    parameter real T_CARRY = 7.8
) (
    input  wire clk,
    input  wire a,
    input  wire b,
    output reg  s,
    output reg  c
);
    reg val_a, val_b;

    initial begin
        val_a = 0;
        val_b = 0;
        s = 0;
        c = 0;
    end

    always @(posedge a) val_a <= 1'b1;
    always @(posedge b) val_b <= 1'b1;

    always @(posedge clk) begin
        if (val_a ^ val_b) begin
            #(T_SUM) s <= 1'b1;
            #2.0 s <= 1'b0;
        end else if (val_a & val_b) begin
            #(T_CARRY) c <= 1'b1;
            #2.0 c <= 1'b0;
        end
        val_a <= 1'b0;
        val_b <= 1'b0;
    end
endmodule
