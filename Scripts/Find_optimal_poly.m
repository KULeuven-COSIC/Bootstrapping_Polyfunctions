load "Scripts/Dependencies.m";

// Legendre expression: the maximum value of e such that p^e divides Factorial(s)
function Legendre(p, s)
    result := 0;
    for i := 1 to Ceiling(Log(p, MaximumOrOne(s))) do
        result +:= Floor(s / (p^i));
    end for;
    return result;
end function;

// Mu function: smallest positive s such that p^e divides Factorial(s)
function Mu(p, e)
    for s := 0 to p^e do
        if Legendre(p, s) ge e then
            return s;
        end if;
    end for;
end function;

// Reduced Mu function: similar function to normal one, but assumes that
// the input is already reduced modulo p^a
// Note: this function does not coincide with the normal Mu function
// for a = 1, because we will later start from a non-linear polynomial
// --> Then it still works because we have Legendre(p, p*s) = s + Legendre(p, s)
function ReducedMu(p, a, e)
    for s := 0 to p^e do
        if a * s + Legendre(p, s) ge e then
            return s;
        end if;
    end for;
end function;

// Return the lowest degree polynomial, with leading coefficient equal
// to p ^ val, that evaluates to zero modulo p ^ e
// The optional parameter a indicates that the input is already reduced
// modulo p^a
function NullPolynomial(p, e: a := 1, val := 0)
    if a eq 1 then
        return (p^val) * (&*[Zx | x - k : k in [0..Mu(p, e - val) - 1]]);
    else
        // Note: it doesn't matter whether we subtract p^a from the polynomial itself or from its variable
        factor := &*[x - k : k in [-((p-1) div 2)..p div 2]];
        return (p^val) * (&*[Zx | factor - k * (p^a) : k in [0..ReducedMu(p, a, e - val) - 1]]);
    end if;
end function;

// Return a matrix whose rows span the null lattice
// The optional parameter a indicates that the input is already reduced
// modulo p^a
function NullLattice(p, e, size: a := 1)
    polynomials := [Zx | ];
    for i := 0 to size - 1 do
        if a eq 1 then
            polynomial := (p ^ (Maximum(0, e - Legendre(p, i)))) * (&*[Zx | x - k : k in [0..i - 1]]);
        else
            // Adapt polynomial in the same way as above
            factor := &*[x - k : k in [-((p-1) div 2)..p div 2]];
            polynomial := (p ^ (Maximum(0, e - a * (i div p) - Legendre(p, i div p)))) *
                          (&*[Zx | factor - k * (p^a) : k in [0..(i div p) - 1]]) * (x ^ (i mod p));
        end if;
        Append(~polynomials, polynomial);
    end for;

    M := ZeroMatrix(Z, size, size);
    for i := 1 to size do
        for j := 1 to size do
            M[i][j] := Coefficient(polynomials[i], j - 1);
        end for;
    end for;
    return M;
end function;



// Given a polynomial defined modulo p^e, return a new polynomial that
// defines the same function but with less terms and smaller coefficients
// - min_exp is the minimum exponent that should occur in the result
// - even_pow indicates if we force the odd exponents to zero
// - odd_pow indicates if we force the even exponents to zero
// - zero_indices determines which other terms are forced to zero
// - lattice indicates if the CVP on the null lattice should be computed
// - min_deg indicates if the resulting polynomial should have minimal degree
// If the lattice flag is set to true, we also give the remaining lattice as a second return value
// It is also possible to pass an optional argument e_inner, then the function assumes that the input
// modulo p^e_inner is already reduced in the balanced set modulo p
function FindOptimalPoly(poly, p, e: min_exp := 0, lattice := true, min_deg := true, even_pow := false, odd_pow := false,
                                     zero_indices := [false : i in [0..Degree(poly)]], e_inner := 1)
    if min_deg then
        for val := 0 to e do
            poly := poly mod NullPolynomial(p, e: a := e_inner, val := val);  // Euclidean division by null polynomials
            poly := poly mod (p^e);
        end for;
    end if;
    size := Degree(poly) + 1;

    // Determine which entries will be forced to zero
    zero_indices := zero_indices[1..size];
    for index := 1 to size div 2 do
        zero_indices[2 * index] or:= even_pow;
    end for;
    for index := 1 to (size + 1) div 2 do
        zero_indices[2 * index - 1] or:= odd_pow;
    end for;
    for index := 1 to Minimum(min_exp, size) do
        zero_indices[index] := true;
    end for;

    // Set up variables for partial Gauss elimination
    seq := [];
    for index := 1 to size do
        if zero_indices[index] then
            seq cat:= [i eq index select 1 else 0 : i in [1..size]];
        end if;
    end for;
    ring := lattice select Z else Integers(p ^ e);  // Ring over which to solve the system of linear equations
    E := Transpose(Matrix(ring, #seq div size, size, seq));

    // Check if system of linear equations is consistent
    null_lattice := ChangeRing(NullLattice(p, e, size: a := e_inner), ring);
    poly_vector := Matrix(ring, 1, size, Eltseq(poly));
    if not IsConsistent(null_lattice * E, -poly_vector * E) then
        error "Optimal polynomial does not exist.";
    end if;

    // Solve system of linear equations
    particular, homogeneous := Solution(null_lattice * E, -poly_vector * E);
    null_space := Matrix(ring, [Eltseq(Basis(homogeneous)[i]) : i in [1..#Basis(homogeneous)]]);

    // Update search space
    poly_vector +:= particular * null_lattice;
    null_lattice := null_space * null_lattice;

    // Early return if lattice is not used
    if not lattice then
        return (Zx!Eltseq(poly_vector)) mod (p ^ e);
    end if;

    // Search for closest vector
    null_vector := ClosestVectors(LatticeWithBasis(null_lattice), Vector(poly_vector))[1];
    return Zx!Eltseq(poly_vector) - Zx!Eltseq(null_vector), null_lattice;
end function;