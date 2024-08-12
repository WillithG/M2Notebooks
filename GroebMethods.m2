-- This notebook contains implementations of the algorithms found in Chapter 2 of Cox, Little, O'Shea
-- The methods found here are not optimized
-- The main purpose is for me to practice M2 and to better understand Groebner bases
-- Any suggestions to improve the code here is appreciated

divAlgo = (f, Lf) -> (
    -- 
    -- f is the polynomial that we want to divide
    -- L is a list of polynomials to divide by
    -- monomial ordering given by the defined ambient ring
    numDivs := #Lf;
    qs := new MutableList from numDivs:0;
    r := 0; -- := is how we do local assignment (is scope not taken care of otherwise?)
    p := f;
    while p != 0 do (
        i := 0;
        divOccurred := false;
        while (i < numDivs and not divOccurred) do (
            -- check if LT(p) divides LT(f)
            if (leadTerm p % leadTerm Lf#i == 0) then (
                -- any easy way to build a polynomial using a list of exponents?
                qs#i = qs#i + leadTerm p // leadTerm Lf#i;
                p = p - (leadTerm p // leadTerm Lf#i) * Lf#i;
                divOccurred = true;
            )
            else (
                i = i + 1;
            );
        );
        if (not divOccurred) then (
            r = r + leadTerm p;
            p = p - leadTerm p;
        );
    );
    return toList(toList qs,r)
);

-- takes two polynomials, outputs their S-polynomial
-- TODO: test this
Spoly = (f,g) -> (
    lcmLMfg := lcm(leadMonomial f, leadMonomial g);
    return lcmLMfg // leadTerm(f) * f - lcmLMfg // leadTerm(g) * g;
)

computeGroeb = (gensL) -> (
    -- gensL :: list of polys which generate an ideal I
    -- output :: list of polys which are a G basis for the input ideal
    G' := gensL;
    gensExtended = true;
    while gensExtended do (
        -- compute all the remainders
        -- get all non-zero remainders
        -- add them to G'
        -- question: should we do them one at a time?
        -- how hard is this going to make refactoring the procedure to be more efficient?
        gensExtended = false;
        allRemainders := for pair in subsets(G', 2) list (divAlgo(Spoly(pair#0, pair#1), G'))#1;
        allNonzeroRemainders := select(allRemainders, i -> i != 0);
        if (#allNonzeroRemainders > 0) then (
            gensExtended = true;
            G' = G' | allNonzeroRemainders;
        ); 
    );
    return G';
)