% A script to import the Recon2.03 file, save it as a .csv for use in
% Python and R.

clc
clear all

load ../data/recon/Recon2.v03.mat
m = modelRecon2beta121114_fixed_updatedID;

% Split the stoichoimetric matrix into irreversible reactions
S = m.S;

% Get rid of highly connected metabolites
met2del = [20,21,26,27,30,31,37,39,40,41,43,44,46,61,64,73,75,77,78,79,80,81,85,86,...
    89,105,130,145,165,201,249,251,259,260,288,295,317,390,392,398,399,412,613,615,839,859,861];
S(met2del,:)=0;

% Make a histogram of the highly connected metabolites
metsums = sum(abs(sign(m.S)),2);
metsums(met2del) = 0;
%hist(metsums,50)
rsums = sum(abs(sign(m.S)),1);

% Get rid of highly connected reactions
r2del = find(rsums > 4 );

r2g = m.rxnGeneMat;
for i = 1:size(m.S,2)
    
    if m.rev(i) == 1
        S(:,size(S,2)+1) = -S(:,i); % Add a row to S
        r2g(size(r2g,1)+1,:) = m.rxnGeneMat(i,:);
    end
    
end

% Save the important files
%csvwrite('../data/recon/S.csv',full(m.S))
cell2csv('../data/recon/genes.csv',m.genes)
%cell2csv('../data/recon/metKeggID.csv',m.metKeggID)
%cell2csv('../data/recon/metCHEBIID.csv',m.metCHEBIID)
%cell2csv('../data/recon/metPubChemID.csv',m.metPubChemID)
%cell2csv('../data/recon/metHMDB.csv',m.metHMDB)
cell2csv('../data/recon/metnames.tsv',m.metNames,'separator', '\t')

% Create a metabolite to gene adjacency amtrix
signmet2gene = sign(S*r2g);
%csvwrite('../data/recon/metGeneMat_sign.csv',full(signmet2gene));

% Repeat but using unsigned data, with highly connected metabolites trimmed
S2 = m.S;
S2(met2del,:) = 0;
S2( : ,r2del) = 0;
unsignmet2genesparse = sign(abs(S2))*m.rxnGeneMat;
%csvwrite('../data/recon/metGeneMat_unsign.csv',full(unsignmet2gene));

g2g = sign(unsignmet2genesparse'*unsignmet2genesparse);
csvwrite('../data/recon/g2g_unsign_sparse.csv',full(g2g));
gsums = sum(g2g,1);
hist(gsums,50)
find(gsums == 0)

S3 = m.S;
S3(met2del,:) = 0;
unsignmet2gene = sign(abs(S3))*m.rxnGeneMat;
%csvwrite('../data/recon/metGeneMat_unsign.csv',full(unsignmet2gene));

g2g = sign(unsignmet2gene'*unsignmet2gene);
csvwrite('../data/recon/g2g_unsign.csv',full(g2g));
gsums2 = sum(g2g,1);
hist(gsums2,50)
find(gsums == 0)
