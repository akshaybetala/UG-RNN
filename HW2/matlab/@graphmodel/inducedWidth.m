function mx = inducedWidth(gm,pi)
% w=inducedWidth(gm,order) : compute the induced width of a graph for a given elimination order
%
 pi = uint32(pi);
 Nfactors = length(gm.factors);
 Nvars = length(pi);
 bucket = cell(1,Nvars); for i=1:Nvars, bucket{i}={}; end;
 [tmp,pir] = sort(pi); pir=uint32(pir);  % compute positions of each var in ordering
 for j=1:Nfactors,
   vtmp = variables(gm.factors{j});
   if (~isempty(vtmp))
     vtmpPos = pir(vtmp);
     vid = min(vtmpPos);
     bucket{vid}{end+1} = vtmp;
   end;
 end;
 mx=0;
 for i=1:Nvars,              % for each variable to eliminate in order
   %fprintf('===Eliminating %d===\n',pi(i));
   vtmp = uint32([]);          %   compute the size of the factor created by elimination:
   for j=1:length(bucket{i}),      %     union all their variables & count vars.
   vtmp=vunion(vtmp,bucket{i}{j});
   end;
   mx = max(mx, length(vtmp)-1);      % 
   vtmp = vdiff(vtmp,pi(i));        % remove the eliminated var
   %fprintf('Leaving '); fprintf('%d ',vtmp); fprintf('\n');
   vtmpPos = pir(vtmp);          % and find the next earliest eliminated var in the bunch
   vid = min(vtmpPos);
   if (vid) bucket{vid}{end+1} = vtmp; end;    % put it in that bucket
 end;

