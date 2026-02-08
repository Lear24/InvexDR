function a = util_evalFnorm(A,B)

PA = util_projection(A);
PB = util_projection(B);

a = norm(PA-PB,'fro');

end 