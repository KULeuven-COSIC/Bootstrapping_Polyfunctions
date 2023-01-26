// Polynomials over integers
Z := Integers();
Zx<x> := PolynomialRing(Z);

// Real and complex numbers
R := RealField(10);
C := ComplexField(10);



// Return the number that is the lowest in absolute value
// If both absolute values are equal, the second one is returned
function MinAbs(a, b)
    if Abs(a) lt Abs(b) then
        return a;
    else
        return b;
    end if;
end function;

// Centered reduction of a sequence mod qi
function CenteredReductionSequence(seq, qi)
    return [Z | MinAbs(coeff mod qi, (coeff mod qi) - qi) : coeff in seq];
end function;

// Centered reduction of a polynomial mod qi
function CenteredReduction(poly, qi)
    return Zx!CenteredReductionSequence(Eltseq(poly), qi);
end function;

// Concatenate zeros to the given array until it reaches the indicated length
function CatZeros(seq, length)
    return seq cat [Z | 0 : i in [1..length - #seq]];
end function;

// Return the maximum value of the sequence (or value) with a minimum of 1
function MaximumOrOne(seq)
    if (Category(seq) eq Category(0)) or (Category(seq) eq Category(R!0)) then
        seq := [seq];
    end if;

    // Compute maximum of extended sequence
    return Maximum(seq cat [1]);
end function;



// Return the lowest digit removal polynomial with respect to the parameters p and e
function GetLowestDigitRemovalPolynomial(p, e)
    // Integers mod p ^ e
    Zpe := Integers(p ^ e);
    Zpe_poly := PolynomialRing(Zpe);
    Zpe_y<y> := quo<Zpe_poly | Zpe_poly.1 ^ ((e - 1) * (p - 1) + 2)>;

    // Compute B(y) and Taylor expansion
    B := ((1 + y) ^ p) - (y ^ p) - 1;
    taylor := Eltseq(p * ((1 + y) ^ p) * (&+[(-B) ^ i : i in [0..(e - 1) * (p - 1) + 1]]));
    taylor := CatZeros(taylor, (e - 1) * (p - 1) + 2);

    // Compute the polynomial
    poly := Zx!0;
    ChangeUniverse(~taylor, Z);
    for m := p to (e - 1) * (p - 1) + 1 do
        factor := taylor[m + 1] / Factorial(m);
        poly +:= &*[x - i : i in [0..m - 1]] * Numerator(factor) * Modinv(Denominator(factor), p ^ e);
    end for;
    return poly mod (p ^ e);
end function;

// Return the lowest digit retain polynomial with respect to the parameters p and e
function GetLowestDigitRetainPolynomial(p, e)
    return (x - GetLowestDigitRemovalPolynomial(p, e)) mod (p ^ e);
end function;