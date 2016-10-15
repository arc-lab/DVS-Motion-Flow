function [vx, vy, It] = computeFlow(x, y, t, pol, N, TH1, TH2, NCOLS, NROWS)

    It_pos=zeros(NROWS,NCOLS);
    It_neg=zeros(NROWS,NCOLS);
    vx=zeros(NROWS,NCOLS); vy=zeros(NROWS,NCOLS);

    % N=5; TH1 = 0.2; TH2 = 0.1; % for real-world seqs

    for ii=1:1:length(t) 
       ptx=x(ii)+1;
       pty=y(ii)+1;
       
%        if ii == 300001
%            keyboard
%        end

       if pol(ii)==1
            It_pos(pty,ptx)=t(ii);
            m=It_pos(max(pty-N,1):min(pty+N, NROWS),max(ptx-N,1):min(ptx+N, NCOLS));
       else
            It_neg(pty,ptx)=t(ii);
            m=It_neg(max(pty-N,1):min(pty+N, NROWS),max(ptx-N,1):min(ptx+N, NCOLS));
       end

       if numel(m) == (2*N+1)*(2*N+1)
           m(abs(m(N+1,N+1)-m)/m(N+1,N+1)>TH1)=0;

           if (sum(m(:)>0))
               [vvx,vvy]=fitplane(m, TH2);
               if isnan(vvx) || isinf(vvx), vvx = 0; end;
               if isnan(vvy) || isinf(vvy), vvy = 0; end;
               
               %vvx(isnan(vvx))= 0; vvy(isnan(vvy)) = 0;
               %vvx(isinf(vvx)) = 0; vvy(isinf(vvy)) = 0;

               if (norm([vvx,vvy])>0)
%                     aa=[vvx(N+1,N+1) vvy(N+1,N+1)];
                    aa=[vvx vvy];
    % %                 vx(pty,ptx)=aa(1)/norm(aa);
    % %                 vy(pty,ptx)=aa(2)/norm(aa);
                    % TODO: x and y are switched somehow
%                     vy(pty,ptx)=aa(1)/norm(aa);
%                     vx(pty,ptx)=aa(2)/norm(aa);

                    vy(pty,ptx)=aa(1);
                    vx(pty,ptx)=aa(2);
               end
           end
       end
    end
% keyboard
   % Return only one time map
   It = max(cat(3, It_pos, It_neg), [], 3);
end
