load "Scripts/Find_optimal_poly.m";

// Test polynomial optimization algorithm on plaintext space mod 2^r
poly := GetLowestDigitRetainPolynomial(2, 16);
new_poly := FindOptimalPoly(poly, 2, 16);
result := true;
for element := 1 to 2^16 do
    if Evaluate(poly, element) mod (2^16) ne Evaluate(new_poly, element) mod (2^16) then
        result := false;
    end if;
end for;
"Test polynomial optimization for plaintext space mod 2^16", result;

// Test minimal power flag
new_poly := FindOptimalPoly(poly, 2, 16: min_exp := 8);
result := true;
for index := 0 to 7 do
    if Coefficient(new_poly, index) ne 0 then
        result := false;
    end if;
end for;
"Test min_exp flag", result;

// Test even power flag
new_poly := FindOptimalPoly(poly, 2, 16: even_pow := true);
result := true;
for index := 1 to 8 do
    if Coefficient(new_poly, 2 * index - 1) ne 0 then
        result := false;
    end if;
end for;
"Test even_pow flag", result;

// Test polynomial optimization algorithm on plaintext space mod 5^r
poly := GetLowestDigitRetainPolynomial(5, 3);
new_poly := FindOptimalPoly(poly, 5, 3);
result := true;
for element := 1 to 5^3 do
    if Evaluate(poly, element) mod (5^3) ne Evaluate(new_poly, element) mod (5^3) then
        result := false;
    end if;
end for;
"Test polynomial optimization for plaintext space mod 5^3", result;

// Test polynomial optimization algorithm for non-minimal degree
poly := GetLowestDigitRetainPolynomial(2, 64);
new_poly := FindOptimalPoly(poly, 2, 35: lattice := false, min_deg := false,
                                         zero_indices := [index mod 8 eq 0 select false else true : index in [0..64]]);



// Test lowest digit retain polynomial for balanced digits instead of non-negative ones
// --> Result is of the form x*f(x^2)
poly := GetLowestDigitRetainPolynomial(17, 4);
poly := (Evaluate(poly, x + (17 div 2)) - (17 div 2)) mod (17^4);
new_poly := FindOptimalPoly(poly, 17, 4: lattice := false);

// Note that the test below does not work for non-negative digit extraction
// --> Working with balanced vs unbalanced digits is not equivalent
// --> Result is of the form x*f(x^4)
poly := x^27;
new_poly := FindOptimalPoly(poly, 3, 4: lattice := false, min_deg := false, zero_indices := [true, true, true,
                                                                                             true, true, false,
                                                                                             true, true, true,
                                                                                             false, true, true,
                                                                                             true, true, true,
                                                                                             true, true, true,
                                                                                             true, true, true,
                                                                                             true, true, true,
                                                                                             true, true, true, true]);



// Test function composition for even p
inner_poly := GetLowestDigitRetainPolynomial(2, 16);
inner_poly := FindOptimalPoly(inner_poly, 2, 16: even_pow := true);
outer_poly := GetLowestDigitRetainPolynomial(2, 64);
outer_poly := FindOptimalPoly(outer_poly, 2, 64: e_inner := 16);

// Check if result is correct
result := true;
for element := 1 to 2^8 do
    bit := Random([0, 1]);
    if Evaluate(outer_poly, element * (2^16) + bit) mod (2^64) ne bit then
        result := false;
    end if;
end for;
"Test outer polynomial of function composition for even p", result;

// Test function composition for odd p
inner_poly := GetLowestDigitRetainPolynomial(3, 8);
inner_poly := (Evaluate(inner_poly, x + (3 div 2)) - (3 div 2)) mod (3^8);
inner_poly := FindOptimalPoly(inner_poly, 3, 8: odd_pow := true);
outer_poly := GetLowestDigitRetainPolynomial(3, 16);
outer_poly := (Evaluate(outer_poly, x + (3 div 2)) - (3 div 2)) mod (3^16);
outer_poly := FindOptimalPoly(outer_poly, 3, 16: odd_pow := true, e_inner := 8);

// Check if result is correct
result := true;
for element := 1 to 2^8 do
    bit := Random([-1, 0, 1]);
    if Evaluate(outer_poly, element * (3^8) + bit) mod (3^16) ne (bit mod (3^16)) then
        result := false;
    end if;
end for;
"Test outer polynomial of function composition for odd p", result;