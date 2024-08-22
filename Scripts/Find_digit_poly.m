load "Scripts/Find_optimal_poly.m";

SetOutputFile("HElib/HElib/src/Magma_define.h": Overwrite := true);

// Generate header file
primes := [ 2, 3, 5, 17, 127 ];
print "long nb_primes = " cat IntegerToString(#primes) cat ";";
print "std::vector<long> primes_list {" cat
                                       &cat[index eq #primes select IntegerToString(primes[index]) else
                                                                    IntegerToString(primes[index]) cat ", " :
                                            index in [1..#primes]] cat "};";

e_inner_list := [ 1, 5, 6, 16 ];
print "long nb_e_inner = " cat IntegerToString(#e_inner_list) cat ";";
print "std::vector<long> e_inner_list {" cat
                                        &cat[index eq #e_inner_list select IntegerToString(e_inner_list[index]) else
                                                                           IntegerToString(e_inner_list[index]) cat ", " :
                                             index in [1..#e_inner_list]] cat "};";

print "std::vector<NTL::ZZX> polynomial_vector[" cat IntegerToString(#primes) cat "][" cat IntegerToString(#e_inner_list) cat "];";

UnsetOutputFile();



SetOutputFile("HElib/HElib/src/Magma_declare.h": Overwrite := true);
print "extern std::vector<NTL::ZZX> polynomial_vector[][" cat IntegerToString(#e_inner_list) cat "];";
UnsetOutputFile();



// Generate all the polynomials
// Use function composition approach: file format is poly<p>_<e_inner>_<e>
for p in primes do
    for e_inner in e_inner_list do
        if ((p eq 17) or (p eq 127)) and (e_inner gt 1) then
            continue;
        end if;
        for e := e_inner + 1 to Round(64 / Log(2, p)) do
            new_e := (p eq 2) and (e mod 2 eq 1) select e + 1 else e;
            poly := GetLowestDigitRetainPolynomial(p, new_e);
            if p ne 2 then
                poly := (Evaluate(poly, x + (p div 2)) - (p div 2)) mod (p^e);
            end if;

            // Make a difference between original version and function composition
            if e_inner eq 1 then
                new_poly := Eltseq(CenteredReduction(FindOptimalPoly(poly, p, e: lattice := (p - 1) * (e - 1) + 1 le 32,
                                                                                 min_deg := p ne 2, even_pow := p eq 2,
                                                                                 odd_pow := p ne 2), p ^ e));
            else
                new_poly := Eltseq(CenteredReduction(FindOptimalPoly(poly, p, e: lattice :=
                                                                   ((p eq 2) and p * ReducedMu(p, e_inner, e) le 10) or
                                                                   ((p eq 3) and p * ReducedMu(p, e_inner, e) le 8) or
                                                                   ((p eq 5) and p * ReducedMu(p, e_inner, e) le 8),
                                                                     odd_pow := p ne 2, e_inner := e_inner), p ^ e));
            end if;

            // Write polynomial to file
            SetOutputFile("HElib/Polynomials/polynomials/poly" cat IntegerToString(p) cat "_" cat IntegerToString(e_inner) cat
                                                                                          "_" cat IntegerToString(e) cat ".txt":
                          Overwrite := true);
            print &cat[index eq #new_poly select IntegerToString(new_poly[index]) else
                                                 IntegerToString(new_poly[index]) cat " " : index in [1..#new_poly]];
            UnsetOutputFile();
        end for;
    end for;
end for;