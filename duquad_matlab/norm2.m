function x = norm2(vec)
x = 0;
for i=1:length(vec)
   x = x + (vec(i)^2); 
end
x = sqrt(x);
end