// This file checks if there exists a bit extraction polynomial that can be
// evaluated recursively with polynomials of lower degree
load "Scripts/Find_optimal_poly.m";

// Check if the given polynomial can be evaluated recursively
function IsRecursivePolynomial(poly)
    // Extract right-hand side of equations
    poly_seq := Eltseq(poly);
    v1 := poly_seq[9]; v2 := poly_seq[7]; v3 := poly_seq[5]; v4 := poly_seq[3];

    // Loop over each possible value for c
    ring := Integers(2^8);
    for c in ring do
        if not IsUnit(c) then
            continue;
        end if;

        a := (ring!v1) / c^2;
        for multiplier := 0 to 1 do
            d := (ring!((v2 / 2) + multiplier * 2^7)) / a / c;
            b := (ring!v3 - a * d^2) / c;

            // Check if result is valid
            if b * d eq ring!v4 then
                return true, ((Z!a)*x^2 + (Z!b)*x), ((Z!c)*x^4 + (Z!d)*x^2);
            end if;
        end for;
    end for;
    return false;
end function;



// Generate all possible polynomials for e = 4
poly := GetLowestDigitRetainPolynomial(2, 4);
new_poly, lattice := FindOptimalPoly(poly, 2, 4: even_pow := true);
set4 := {};
for it := 1 to 10000 do
    Include(~set4, new_poly);

    // Generate new random polynomial using the lattice
    vector := [Random(2^4 - 1) : i in [1..NumberOfRows(lattice)]];
    vector := Matrix(Z, 1, #vector, vector);
    new_poly := (new_poly + Zx!Eltseq(vector * lattice)) mod 2^4;
end for;



// Generate all possible polynomials for e = 8 and check if they are easy to evaluate
poly := GetLowestDigitRetainPolynomial(2, 8);
new_poly, lattice := FindOptimalPoly(poly, 2, 8: even_pow := true);
set8 := {};
for it := 1 to 10000 do
    if IsRecursivePolynomial(new_poly) then
        _, f1, f2 := IsRecursivePolynomial(new_poly);
        Include(~set8, <f1, f2, Evaluate(f1, f2)>);
    end if;

    // Generate new random polynomial using the lattice
    vector := [Random(2^8 - 1) : i in [1..NumberOfRows(lattice)]];
    vector := Matrix(Z, 1, #vector, vector);
    new_poly := (new_poly + Zx!Eltseq(vector * lattice)) mod 2^8;
end for;