load "Scripts/Find_optimal_poly.m";

// Regularize the given set of univariate polynomials by dividing by the
// greatest common divisor of all coefficients
function RegularizePolynomials(polynomials)
    // Convert data structures
    _, p, e := IsPrimePower(#BaseRing(polynomials[1]));
    ChangeUniverse(~polynomials, Zx);

    // Find maximum valuation and use it for division
    val := Valuation(GCD(&cat[Eltseq(polynomial) : polynomial in polynomials]), p);
    return [PolynomialRing(Integers(p ^ (e - val))) | polynomial div (p ^ val) : polynomial in polynomials];
end function;

// Find one common root of the given univariate polynomials modulo p ^ e
// If no root can be found after nb_iteration attempts, we return -1
function FindRoot(polynomials: nb_iterations := 1000)
    // Divide by greatest common divisor of all coefficients
    polynomials := RegularizePolynomials(polynomials);

    // Construct all necessary structures
    _, p, e := IsPrimePower(#BaseRing(polynomials[1]));
    Fp := GF(p);
    Fp_y<y> := PolynomialRing(Fp);

    // First find the roots modulo p: ideal is principal over field
    poly := GCD([Evaluate(Zx!polynomial, y) : polynomial in polynomials]);
    roots := (poly eq 0) select {Z | element : element in Fp} else {Z | element[1] : element in Roots(poly)};

    // Return -1 if no root exists
    if #roots eq 0 then
        return -1;
    end if;

    // Perform a number of random attempts to Hensel lift the result
    for iteration := 1 to nb_iterations do
        root := Random(roots); prec := 1;
        while prec lt e do
            // Compute new polynomial equations for random root
            extra_prec := Minimum(e - prec, prec);
            new_polynomials := [(Evaluate(Zx!polynomial, (p ^ prec) * x + root) div (p ^ prec)) mod (p ^ extra_prec) :
                                polynomial in polynomials];
            new_polynomials := [CatZeros(Eltseq(polynomial), 2) : polynomial in new_polynomials];

            // System of equations
            LHS := [new_polynomials[index][2] : index in [1..#new_polynomials]];
            RHS := [-new_polynomials[index][1] : index in [1..#new_polynomials]];
            if not IsConsistent(Matrix(Integers(p ^ extra_prec), 1, #new_polynomials, LHS),
                                Matrix(Integers(p ^ extra_prec), 1, #new_polynomials, RHS)) then
                break;
            end if;

            // Find solutions: homogeneous part can only be of dimension 1
            result, homogeneous := Solution(Matrix(Integers(p ^ extra_prec), 1, #new_polynomials, LHS),
                                            Matrix(Integers(p ^ extra_prec), 1, #new_polynomials, RHS));
            homogeneous := Basis(homogeneous); homogeneous := homogeneous eq [] select 0 else Z!homogeneous[1][1];
            result := Z!Eltseq(result)[1];

            // Choose one random solution to continue the search
            root := ((p ^ prec) * (result + Random(p ^ extra_prec - 1) * homogeneous) + root) mod (p ^ (prec + extra_prec));
            prec +:= extra_prec;
        end while;

        // Result is only valid if it has sufficiently high precision
        if prec ge e then
            return root;
        end if;
    end for;
    return -1;
end function;

// Return the minimum value of the sequence with a maximum of the given default
function MinimumOrDefault(seq, defaultValue)
    return Minimum(seq cat [defaultValue]);
end function;

// Convert a given element of a multivariate polynomial ring that is in fact
// univariate into a true univariate polynomial
function GetUnivariatePolynomial(poly, currentVariable)
    _, poly := IsUnivariate(poly, currentVariable);
    return poly;
end function;

// Perform back substitution on the given Gröbner basis
// This function returns at most one random solution
function BackSubstitutionBody(G, N: currentVariable := N)
    // Solution found
    if #G eq 0 then
        return [[0 : i in [1..currentVariable]]];
    elif currentVariable eq 0 then
        for element in G do
            if element ne 0 then
                return [];
            end if;
        end for;
        return [[]];
    end if;

    // Start with last equation
    if G[#G] eq 0 then
        return BackSubstitutionBody(Remove(G, #G), N: currentVariable := currentVariable);
    end if;

    // Compute result via recursive call
    FpN := Parent(G[1]);
    firstIndex := MinimumOrDefault([index : index in [1..#G] | &and[IsUnivariate(G[next_index], currentVariable) :
                                                                    next_index in [index..#G]]], #G);
    root := IsUnivariate(G[firstIndex], currentVariable) select
            FindRoot([GetUnivariatePolynomial(G[index], currentVariable) : index in [firstIndex..#G]]) else Z!Random(BaseRing(FpN));
    
    // Return empty sequence if there is no solution
    if root eq -1 then
        return [];
    end if;

    // Compute result recursively assuming that the given root is correct
    recursive := BackSubstitutionBody([Evaluate(equation, [FpN | index eq currentVariable select root else FpN.index :
                                                           index in [1..N]]) : equation in G], N:
                                                           currentVariable := currentVariable - 1);
    
    // Extend recursive result with current root
    if recursive ne [] then
        return [Append(recursive[1], root)];
    end if;
    return [];
end function;

// Perform back substitution on the given Gröbner basis
// This function iterates until a valid solution is found
function BackSubstitution(G, N)
    while true do
        result := BackSubstitutionBody(G, N: currentVariable := N);
        if result ne [] then
            return result[1];
        end if;
    end while;
end function;

// Get the ideal corresponding to the given polynomial optimization problem
// The ideal is computed over the integers modulo p ^ e
function GetIdeal(poly, lattice, previousPoly, p, e)
    // Construct all necessary structures
    N := 3 * #previousPoly + NumberOfRows(lattice);
    Z_poly := PolynomialRing(Integers(p ^ e), N);
    Z_poly_y<y> := PolynomialRing(Z_poly);

    // Construct equations
    LHS := Vector(Eltseq(poly)) + Vector([Z_poly.index : index in [3 * #previousPoly + 1..N]]) *
           Matrix(Z_poly, NumberOfRows(lattice), NumberOfColumns(lattice), Eltseq(lattice));
    RHS := &+[Z_poly.index * Evaluate(previousPoly[index], y) : index in [1..#previousPoly]] -
          (&+[Z_poly.index * Evaluate(previousPoly[index - #previousPoly], y) :
              index in [#previousPoly + 1..2 * #previousPoly]]) *
          (&+[Z_poly.index * Evaluate(previousPoly[index - 2 * #previousPoly], y) :
              index in [2 * #previousPoly + 1..3 * #previousPoly]]);

    // Convert to regular sequence
    LHS := Eltseq(LHS);
    RHS := Eltseq(RHS);

    // Construct ideal
    return [LHS[index] - RHS[index] : index in [1..#LHS]] cat [Z_poly.1, Z_poly.3, Z_poly.4, Z_poly.9 + 1];
    //return [LHS[index] - RHS[index] : index in [1..#LHS]];    // Normal procedure (if code example below is not used)
end function;

// Perform Gröbner basis analysis for recursive evaluation of the given polynomial
function FindRecursivePolynomial(poly, lattice, previousPoly, p, e)
    N := 3 * #previousPoly + NumberOfRows(lattice);
    I := GetIdeal(poly, lattice, previousPoly, p, e);
    G := GroebnerBasis(I);
    result := BackSubstitution(G, N);

    // Create polynomial based on result
    return (&+[result[index] * previousPoly[index] : index in [1..#previousPoly]] -
           (&+[result[index] * previousPoly[index - #previousPoly] :
              index in [#previousPoly + 1..2 * #previousPoly]]) *
           (&+[result[index] * previousPoly[index - 2 * #previousPoly] :
              index in [2 * #previousPoly + 1..3 * #previousPoly]])) mod (p ^ e),
           CenteredReductionSequence(result[1..3 * #previousPoly], p ^ e);
end function;



// Find optimal polynomial using Gröbner basis analysis
poly := GetLowestDigitRetainPolynomial(2, 16);
new_poly, lattice := FindOptimalPoly(poly, 2, 16: even_pow := true);
rec_poly, result := FindRecursivePolynomial(new_poly, lattice, [x^2, x^4, 14641*x^8 + 22748*x^6 + 8836*x^4 + 112*x^2], 2, 16);

// Recursive evaluation of lowest digit retain polynomial (as implemented in HElib):
// * f2 = x^2
// * f4 = (f2)^2
// * f8 = 112*f2 + (94*f2 + 121*f4)^2
// * f16 = 11136*f4 - (15364*f4 - 14115*f8) * (28504*f2 + 8968*f4 - f8)
// ==> Note that reduction is allowed modulo 2^16 only

// Concrete values (as implemented in HElib):
// * f2 = x^2
// * f4 = x^4
// * f8 = 14641*x^8 + 22748*x^6 + 8836*x^4 + 112*x^2
// * f16 = 10941*x^16 + 31832*x^14 - 28268*x^12 - 13992*x^10 - 19072*x^8 + 23680*x^6 - 5120*x^4