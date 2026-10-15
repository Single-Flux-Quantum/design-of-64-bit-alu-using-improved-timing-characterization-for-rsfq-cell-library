`timescale 1ps/100fs

// Verification Testbench for 64-Bit RSFQ Wave-Pipelined ALU
// Reference: Inamdar et al., IEEE Trans. Appl. Supercond., Vol. 31, No. 5, August 2021

module tb_alu64;
    reg         clk;
    reg  [2:0]  opcode;
    reg  [63:0] a;
    reg  [63:0] b;
    reg         cin;
    wire [63:0] result;
    wire        cout;

    reg  [63:0] expected_res;
    reg         expected_cout;

    integer pass_count;
    integer fail_count;

    // Instantiate 64-Bit ALU
    alu64_top u_alu (
        .clk(clk),
        .opcode(opcode),
        .a(a),
        .b(b),
        .cin(cin),
        .result(result),
        .cout(cout)
    );

    // Clock Generator: 25 GHz Nominal (Period = 40 ps)
    always #20.0 clk = ~clk;

    initial begin
        clk = 0;
        opcode = 0;
        a = 0;
        b = 0;
        cin = 0;
        pass_count = 0;
        fail_count = 0;

        #100;
        $display("================================================================================");
        $display("TESTBENCH: 64-BIT RSFQ ARITHMETIC LOGIC UNIT (INAMDAR ET AL. 2021)");
        $display("================================================================================");

        // ---------------------------------------------------------------------
        // Test 1: 64-bit ADD Operation
        // ---------------------------------------------------------------------
        @(posedge clk);
        opcode = 3'b000;
        a = 64'h12345678_9ABCDEF0;
        b = 64'h0FEDCBA9_87654321;
        cin = 1'b0;
        expected_res = a + b + cin;
        @(posedge clk);
        @(posedge clk);
        #10;
        if (result === expected_res) begin
            $display("[PASS] 64-Bit ADD: 0x%016h + 0x%016h = 0x%016h", a, b, result);
            pass_count = pass_count + 1;
        end else begin
            $display("[FAIL] 64-Bit ADD: Expected 0x%016h, got 0x%016h", expected_res, result);
            fail_count = fail_count + 1;
        end

        // ---------------------------------------------------------------------
        // Test 2: 64-bit SUB Operation
        // ---------------------------------------------------------------------
        @(posedge clk);
        opcode = 3'b001;
        a = 64'hFEDCBA98_76543210;
        b = 64'h12345678_90ABCDEF;
        cin = 1'b1;
        expected_res = a - b;
        @(posedge clk);
        @(posedge clk);
        #10;
        if (result === expected_res) begin
            $display("[PASS] 64-Bit SUB: 0x%016h - 0x%016h = 0x%016h", a, b, result);
            pass_count = pass_count + 1;
        end else begin
            $display("[FAIL] 64-Bit SUB: Expected 0x%016h, got 0x%016h", expected_res, result);
            fail_count = fail_count + 1;
        end

        // ---------------------------------------------------------------------
        // Test 3: 64-bit bitwise AND
        // ---------------------------------------------------------------------
        @(posedge clk);
        opcode = 3'b010;
        a = 64'hF0F0F0F0_AAAAAAAA;
        b = 64'hFFFF0000_55555555;
        expected_res = a & b;
        @(posedge clk);
        @(posedge clk);
        #10;
        if (result === expected_res) begin
            $display("[PASS] 64-Bit AND: 0x%016h & 0x%016h = 0x%016h", a, b, result);
            pass_count = pass_count + 1;
        end else begin
            $display("[FAIL] 64-Bit AND: Expected 0x%016h, got 0x%016h", expected_res, result);
            fail_count = fail_count + 1;
        end

        // ---------------------------------------------------------------------
        // Test 4: 64-bit bitwise OR
        // ---------------------------------------------------------------------
        @(posedge clk);
        opcode = 3'b011;
        a = 64'hF0F0F0F0_AAAAAAAA;
        b = 64'h0F0F0F0F_55555555;
        expected_res = a | b;
        @(posedge clk);
        @(posedge clk);
        #10;
        if (result === expected_res) begin
            $display("[PASS] 64-Bit OR: 0x%016h | 0x%016h = 0x%016h", a, b, result);
            pass_count = pass_count + 1;
        end else begin
            $display("[FAIL] 64-Bit OR: Expected 0x%016h, got 0x%016h", expected_res, result);
            fail_count = fail_count + 1;
        end

        // ---------------------------------------------------------------------
        // Test 5: 64-bit bitwise XOR
        // ---------------------------------------------------------------------
        @(posedge clk);
        opcode = 3'b100;
        a = 64'hFFFFFFFF_00000000;
        b = 64'hAAAA5555_AAAA5555;
        expected_res = a ^ b;
        @(posedge clk);
        @(posedge clk);
        #10;
        if (result === expected_res) begin
            $display("[PASS] 64-Bit XOR: 0x%016h ^ 0x%016h = 0x%016h", a, b, result);
            pass_count = pass_count + 1;
        end else begin
            $display("[FAIL] 64-Bit XOR: Expected 0x%016h, got 0x%016h", expected_res, result);
            fail_count = fail_count + 1;
        end

        $display("================================================================================");
        $display("TESTBENCH SUMMARY: %0d Passed, %0d Failed", pass_count, fail_count);
        $display("================================================================================");
        #100;
        $finish;
    end
endmodule
