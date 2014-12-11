function fid = latexTab(M, decimals, filename, mode, label)
%
% M: matrix, desimals: int, filename: 'name.txt', mode: eg 'w', label: '%something'.
%
fid = fopen(filename, mode);
[r c] = size(M);
fprintf(fid, '%s: \n',label);
for i = 1:r
    for j = 1:c
        %fprintf(fid, '%.4g ', M(i,j));   % se http://amath.colorado.edu/computing/Matlab/OldTechDocs/ref/fprintf.html
        if M(i,j) > 10^(-4) && M(i,j) < 10^(-2)
            fprintf(fid, strcat('%f '), M(i,j));
        else
            fprintf(fid, strcat('%f '), M(i,j));
        end
        if j ~= c
            fprintf(fid, '& ');
        end
    end
    fprintf(fid, '\\\\ \n');
end
fprintf(fid, '\n');
fclose(fid);
end

