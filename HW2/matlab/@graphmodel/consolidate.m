function gm=consolidate(gm)

% see if each factor i can be joined with (one or more) others
for i=2:length(gm.factors),  %if (~isempty(gm.factors{i})),
  vars=variables(gm.factors{i});
  if (isempty(vars)) 
    gm.factors{1}=gm.factors{1}*gm.factors{i};
  else
    candid = vdiff( gm.vNbrs{vars(1)} ,uint32(i) );                % intersecting variable lists for candidates
    for j=vars(2:end), candid=vintersect(candid,gm.vNbrs{j}); end;
    if (length(candid)>0)                                          % if so, add a portion of factor i
      F = gm.factors{i}.^(1/length(candid));                      %  to each of the others
      for j=candid, gm.factors{j}=gm.factors{j} * F; end;
      gm = removeFactor(gm,i);                                    % then get rid of factor i
    end;
  end; 
end;

% push non-empty factors to the beginning of the list
gm.vacant = sort(gm.vacant,2,'descend');
i=length(gm.factors); j=1; k=gm.nVacant;
while (j<k)
  while (i<=gm.vacant(j)) 
    if (i==gm.vacant(j)) i=i-1; end; 
    j=j+1; 
  end;
  if (i>gm.vacant(k))
    gm.factors{gm.vacant(k)}=gm.factors{i};
    gm.factors{i}=factor();
    gm.vacant(k)=i;
    k=k-1;
    i=i-1;
  end;
end;
gm.vacant = sort(gm.vacant,2,'descend');

%for i=1:length(gm.factors),
%  tmp=gm.factors{i};
%  if (~isempty(tmp))
%    gm=removeFactor(gm,i);
%    gm=addFactor(gm,tmp);
%  end;
%end;
