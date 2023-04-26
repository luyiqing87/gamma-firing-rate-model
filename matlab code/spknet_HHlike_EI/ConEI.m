function ConMat = ConEI(Ne,Ni,pee,pii,pei,pie)
% Generate a connectivity matrix for a network with Ne excitatory cells, Ni
% inhibitory cells and connectivity probability pee, pii, pei, pie
% Example: 
% ConMat=ConEI(800,200,0.05,0.3,0.3,0.3);

MatEE=zeros(Ne);
MatII=zeros(Ni);
MatEI=zeros(Ne,Ni);
MatIE=zeros(Ni,Ne);

RandMat=rand(Ne+Ni);
for i=1:Ne
    for j=1:Ne
        if RandMat(i,j) <= pee
            MatEE(i,j)=1;
        end
    end
end

for i=1:Ne
    for j=1:Ni
        if RandMat(i,j+Ne) <= pei
            MatEI(i,j)=1;
        end
    end
end

for i=1:Ni
    for j=1:Ne
        if RandMat(i+Ne,j) <= pie
            MatIE(i,j)=1;
        end
    end
end

for i=1:Ni
    for j=1:Ni
        if RandMat(i+Ne,j+Ne) <= pii
            MatII(i,j)=1;
        end
    end
end

ConMat=[MatEE,MatEI;MatIE,MatII];
for i=1:(Ne+Ni)
    ConMat(i,i)=0;
end

end
