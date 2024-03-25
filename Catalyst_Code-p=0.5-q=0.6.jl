#This function sorts the vector in descending order
function _vsort(p::Vector)
    return sort(p,rev=true)
end

#This function pads to length d
function _vpad(d :: Int, p::Vector)
    pad = zeros(d)
    pad[1:length(p)] = p
    return pad
end

#This function computes the classical fidelity of two vectors (i.e. Bhattacharya distance squared)
function _cFidelity(p::Vector, q::Vector)
    if length(p) > length(q)
        d = length(p)
        q = _vpad(d,q)
    elseif length(q) < length(p)
        d = length(q)
        p = _vpad(d,p)
    else
        d = length(q)
    end

    F = 0
    for i in 1:d
        F += sqrt(p[i]*q[i])
    end
    return F^2
end

#This returns the catalyst fidelity for a bernoulli catalyst
function _bernFidelity(p,q,r)
    p = [p,1-p]
    q = [q,1-q]
    r = [r,1-r]
    pr=_vsort(kron(p,r))
    qr=_vsort(kron(q,r))
    return _cFidelity(pr,qr)
end

#This generates a random distribution for a fixed dimension
function _randDist(dim)
    distVect = rand(Float64,(dim,1))
    return distVect/sum(distVect)
end

#Hardcode p and q at the start
pval = 0.5
qval = 0.6

println(_bernFidelity(pval,qval,0))


#d = 2
pvec = [pval,1-pval]
qvec = [qval,1-qval]
dim = 2
bestR= zeros(size(dim))
currentF = 0
maxF = 0
ctr = 0
step_size = 0.005
for r1 in [0:step_size:0.5;]
    for r2 in [r1:step_size:0.5;]
        r = [r1,r2]
        if ctr % 100000 == 0
            println("Now on ", r)
        end
        ctr += 1
        if (r1+r2) == 1
            pr=_vsort(kron(pvec,r))
            qr=_vsort(kron(qvec,r))
            currentF = _cFidelity(pr,qr)
            if currentF > maxF
                maxF = currentF
                bestR = r
            else
            end
        end
    end
end



#d=3
pvec = [pval,1-pval]
qvec = [qval,1-qval]
dim = 3
bestR= zeros(size(dim))
currentF = 0
maxF = 0
ctr = 0
step_size = 0.005 
for r1 in [0:step_size:0.4;]
    for r2 in [r1:step_size:0.5;]
        for r3 in [r2:step_size:0.5;]
            r = [r1,r2,r3]
            if ctr % 10^8 == 0
                println("Now on ", r)
                ctr = 0
            end
            ctr += 1
            if (r1+r2+r3) == 1
                pr=_vsort(kron(pvec,r))
                qr=_vsort(kron(qvec,r))
                currentF = _cFidelity(pr,qr)
                if currentF > maxF
                    maxF = currentF
                    bestR = r
                else
                end
            end
        end
    end
end


#d=4
pvec = [pval,1-pval]
qvec = [qval,1-qval]
dim = 4
bestR= zeros(size(dim))
currentF = 0
maxF = 0
ctr = 0
step_size = 0.005 
for r1 in [0:step_size:0.4;]
    for r2 in [r1:step_size:0.4;]
        for r3 in [r2:step_size:0.5;]
            for r4 in [r3:step_size:0.5;]
                r = [r1,r2,r3,r4]
                if ctr % 10^8 == 0
                    println("Now on ", r)
                    ctr = 0
                end
                ctr += 1
                if (r1+r2+r3+r4) == 1
                    pr=_vsort(kron(pvec,r))
                    qr=_vsort(kron(qvec,r))
                    currentF = _cFidelity(pr,qr)
                    if currentF > maxF
                        maxF = currentF
                        bestR = r
                    else
                    end
                end
            end
        end
    end
end

#d=5
pvec = [pval,1-pval]
qvec = [qval,1-qval]
dim = 5
bestR= zeros(size(dim))
currentF = 0
maxF = 0
ctr = 0
step_size = 0.005 
for r1 in [0:step_size:0.3;]
    for r2 in [r1:step_size:0.4;]
        for r3 in [r2:step_size:0.4;]
            for r4 in [r3:step_size:0.5;]
                for r5 in [r4:step_size:0.5;]
                    r = [r1,r2,r3,r4,r5]
                    if ctr % 10^8 == 0
                        println("Now on ", r)
                        ctr = 0
                    end
                    ctr += 1
                    if sum(r) == 1
                        pr=_vsort(kron(pvec,r))
                        qr=_vsort(kron(qvec,r))
                        currentF = _cFidelity(pr,qr)
                        if currentF > maxF
                            maxF = currentF
                            bestR = r
                        else
                        end
                    end
                end
            end
        end
    end
end

#d=6
pvec = [pval,1-pval]
qvec = [qval,1-qval]
dim = 6
bestR= zeros(size(dim))
currentF = 0
maxF = 0
ctr = 0
step_size = 0.005
for r1 in [0:step_size:0.3;]
    for r2 in [r1:step_size:0.3;]
        for r3 in [r2:step_size:0.4;]
            for r4 in [r3:step_size:0.4;]
                for r5 in [r4:step_size:0.5;]
                    for r6 in [r5:step_size:0.5;]
                        r = [r1,r2,r3,r4,r5,r6]
                        if ctr % 10^8 == 0
                            println("Now on ", r)
                            ctr = 0
                        end
                        ctr += 1
                        if sum(r) == 1
                            pr=_vsort(kron(pvec,r))
                            qr=_vsort(kron(qvec,r))
                            currentF = _cFidelity(pr,qr)
                            if currentF > maxF
                                maxF = currentF
                                bestR = r
                            else
                            end
                        end
                    end
                end
            end
        end
    end
end

#d=7
pvec = [pval,1-pval]
qvec = [qval,1-qval]
dim = 7
bestR= zeros(size(dim))
currentF = 0
maxF = 0
ctr = 0
step_size = 0.005
for r1 in [0.0:step_size:0.15;]
    for r2 in [r1:step_size:0.15;]
        for r3 in [r2:step_size:0.17;]
            for r4 in [r3:step_size:0.19;]
                for r5 in [r4:step_size:0.23;]
                    for r6 in [r5:step_size:0.27;]
                        for r7 in [r6:step_size:0.31;]
                            r = [r1, r2, r3, r4, r5, r6, r7]
                            if ctr % 10^8 == 0
                                println("Now on ", r)
                                ctr = 0
                            end
                            ctr += 1
                            if sum(r) == 1
                                pr=_vsort(kron(pvec,r))
                                qr=_vsort(kron(qvec,r))
                                currentF = _cFidelity(pr, qr)
                                if currentF > maxF
                                    maxF = currentF
                                    bestR = r
                                else
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

#d=8
pvec = [pval,1-pval]
qvec = [qval,1-qval]
dim = 8
bestR= zeros(size(dim))
currentF = 0
maxF = 0
ctr = 0
step_size = 0.005
for r1 in [0:step_size:0.13;]
    for r2 in [r1:step_size:0.13;]
        for r3 in [r2:step_size:0.14;]
            for r4 in [r3:step_size:0.16;]
                for r5 in [r4:step_size:0.18;]
                    for r6 in [r5:step_size:0.21;]
                        for r7 in [r6:step_size:0.25;]
                            for r8 in [r7:step_size:0.29;]
                                r = [r1,r2,r3,r4,r5,r6,r7,r8]
                                if ctr % 10^8 == 0
                                    println("Now on ", r)
                                    ctr = 0
                                end
                                ctr += 1
                                if sum(r) == 1
                                    pr=_vsort(kron(pvec,r))
                                    qr=_vsort(kron(qvec,r))
                                    currentF = _cFidelity(pr,qr)
                                    if currentF > maxF
                                        maxF = currentF
                                        bestR = r
                                    else
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
