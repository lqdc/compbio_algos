%% estatic_matrix.m
%% 
%% Linear Algebraic applications to computing long-range electrostatic
%% interactions between atoms.

numCharges = 10;
Nsep = 50;
sepX = logspace(-5,5,Nsep);
separation = zeros(size(sepX,2),3);
separation(:,1) = sepX';
pointCharges = zeros(2 * numCharges, 3);
eps = 2.397e-4;          % epsilon naught in AKMA units
q = rand(numCharges,1);

nTrial = 100;
for iter=1:Nsep
 r1error(iter) =  0;
 r2error(iter) =  0;
 r8error(iter) = 0;

 %for trial=1:nTrial
  %% produce a list of source charge points and destination charge
  %% points  -- be sure to select symmetric (for eigenvalues)  or 
  %% assymetric (for the svd) set.  
  for i=1:numCharges
	 srcCharges(i,:) = rand(1,3);
         % symmetric charge points: 
	 %destCharges(i,:) = reflect(srcCharges(i,:), sepX(iter));
         % assymentic charge points:
	 destCharges(i,:) = rand(1,3) + sepX(iter);          
  end

  
  symmetric_decomposition = 0;
  P = zeros(numCharges,numCharges);
  for i=1:numCharges
      for j=1:numCharges
          if i ~= j
            P(i,j) = 1/(2*eps*4*pi)/norm(srcCharges(i)-destCharges(j));
          end
      end
  end
  if symmetric_decomposition
  %%checking for symmetricity
      for i=1:numCharges
          for j=1:numCharges
              if str2num(sprintf('%5.4f',P(i,j))) ~= str2num(sprintf('%5.4f',P(j,i)))
                  sprintf('at i, j it is %5.2f',P(i,j));
                  sprintf('at j, i it is %5.2f',P(j,i));
                  display(sprintf('not symmetric at %d %d',i,j));
              end
          end
      end
  end
  %% decompose A into eigenvalues and eigenvectors using the
  %% eig command (or into singular vectors, using the svd command)
  if symmetric_decomposition == 1
      [V,D] = eig(P);
      %round matrices and check for equality
      newP = round(V*D*V' .* 10^4) ./ 10^4;
      P = round(P .* 10^4) ./ 10^4;      
      if P ~= newP
          display(sprintf('matrices are different'))
      end

      
      %checking if orthonormal
      not_orthonormal = 0;
      for i=1:size(V,1)
          for j=1:size(V,1)
            dot_product = dot(V(:,i),V(:,j));
            if (i == j) && (round(dot_product) ~= 1)
                not_orthonormal = 1;
                display(sprintf('dot product at i = %d j = %d is %f and should be 1',i,j, dot_product))
            elseif (i ~= j) && (round(dot_product) ~= 0)
                not_orthonormal = 1;
                display(sprintf('dot product at i = %d j = %d is %f and should be 0',i, j, dot_product))

            end
          end
      end
      if not_orthonormal == 0
          %display(sprintf('matrix V is orthonormal'));
      else
          display(sprintf('matrix V is not orthonormal'));
          P
          D
          V
      end
      display(sprintf('separation distance is %f',sepX(iter)))
  else
      [U, S, V] = svd(P);
      newP = round(U*S*V' .* 10^4) ./ 10^4;
      P = round(P .* 10^4) ./ 10^4;
      if P ~= newP
          display(sprintf('matrices are different'))
      end
      %checking if V is orthonormal
      not_orthonormal = 0;
      for i=1:size(V,1)
          for j=1:size(V,1)
            dot_product = dot(V(:,i),V(:,j));
            if (i == j) && (round(dot_product) ~= 1)
                not_orthonormal = 1;
                display(sprintf('dot product at i = %d j = %d is %f and should be 1',i,j, dot_product))
            elseif (i ~= j) && (round(dot_product) ~= 0)
                not_orthonormal = 1;
                display(sprintf('dot product at i = %d j = %d is %f and should be 0',i, j, dot_product))

            end
          end
      end
      %checking if U is orthonormal
      not_orthonormal = 0;
      for i=1:size(U,1)
          for j=1:size(U,1)
            dot_product = dot(U(:,i),U(:,j));
            if (i == j) && (round(dot_product) ~= 1)
                not_orthonormal = 1;
                display(sprintf('dot product at i = %d j = %d is %f and should be 1',i,j, dot_product))
            elseif (i ~= j) && (round(dot_product) ~= 0)
                not_orthonormal = 1;
                display(sprintf('dot product at i = %d j = %d is %f and should be 0',i, j, dot_product))
            end
          end
      end
  end
  
  %% create rank one, two, 8 approximations to A.
 if symmetric_decomposition == 1
    D_sorted=diag(sort(diag(D),'descend')); % make diagonal matrix out of sorted diagonal values of input D
    [c, ind]=sort(diag(D),'descend'); % store the indices of which columns the sorted eigenvalues come from
    V_sorted=V(:,ind); % arrange the columns in this order
    q_sorted=q(ind);%sort the q vector as well
    
    %now make rank approximations of the matrix
    D_k1 = zeros(1, 1);
    D_k2 = zeros(2, 2);
    D_k8 = zeros(8, 8);
    D_k1(1:1,1:1) = D_sorted(1:1,1:1);
    D_k2(1:2,1:2) = D_sorted(1:2,1:2);
    D_k8(1:8,1:8) = D_sorted(1:8,1:8);
    
    V_k1(:,1:1) = V_sorted(:,1:1);
    V_k2(:,1:2) = V_sorted(:,1:2);
    V_k8(:,1:8) = V_sorted(:,1:8);
    
    P_prev = V_sorted * D_sorted * V_sorted';
    P_approx1 = V_k1 * D_k1 * V_k1';
    P_approx2 = V_k2 * D_k2 * V_k2';
    P_approx8 = V_k8 * D_k8 * V_k8';

    
    decomposedP = P_prev * q_sorted;
    approxP8 = P_approx8 * q_sorted;
    approxP2 = P_approx2 * q_sorted;
    approxP1 = P_approx1 * q_sorted;
    realP = P*q;
    
    P_difference1(iter) = norm(decomposedP - approxP1);
    P_difference2(iter) = norm(decomposedP - approxP2);
    P_difference8(iter) = norm(decomposedP - approxP8);
 else
    S_sorted=diag(sort(diag(S),'descend')); % make diagonal matrix out of sorted diagonal values of input D
    [c, ind]=sort(diag(S),'descend'); % store the indices of which columns the sorted eigenvalues come from
    V_sorted=V(:,ind); % arrange the columns in this order
    U_sorted=U(:,ind);
    q_sorted=q(ind);%sort the q vector as well
    
    %now make rank approximations of the matrix
    S_k1(1:1,1:1) = S_sorted(1:1,1:1);
    S_k2(1:2,1:2) = S_sorted(1:2,1:2);
    S_k8(1:8,1:8) = S_sorted(1:8,1:8);
    
    V_k1(:,1:1) = V_sorted(:,1:1);
    V_k2(:,1:2) = V_sorted(:,1:2);
    V_k8(:,1:8) = V_sorted(:,1:8);

    U_k1(:,1:1) = U_sorted(:,1:1);
    U_k2(:,1:2) = U_sorted(:,1:2);
    U_k8(:,1:8) = U_sorted(:,1:8);
    
    P_prev = U_sorted * S_sorted * V_sorted';
    P_approx1 = U_k1 * S_k1 * V_k1';
    P_approx2 = U_k2 * S_k2 * V_k2';
    P_approx8 = U_k8 * S_k8 * V_k8';

    
    decomposedP = P_prev * q_sorted;
    approxP8 = P_approx8 * q_sorted;
    approxP2 = P_approx2 * q_sorted;
    approxP1 = P_approx1 * q_sorted;
    realP = P*q;

    P_difference1(iter) = norm(decomposedP - approxP1);
    P_difference2(iter) = norm(decomposedP - approxP2);
    P_difference8(iter) = norm(decomposedP - approxP8);
 end
  %% see how good the low rank approximation is (using 
  %% the random charge vector q (use variable error1, error2,
  %% and error8)

  %r1error(iter) =  r1error(iter) + error1;
  %r2error(iter) =  r2error(iter) + error2;
  %r8error(iter) =  r8error(iter) + error8;

 %end
  
 %r1error(iter) =  r1error(iter) / nTrial;
 %r2error(iter) =  r2error(iter) / nTrial;
 %r8error(iter) =  r8error(iter) / nTrial;

end
%{
subplot(2,2,1)
title('Cutoff == 1')
hold on
subplot(2,2,3)
title('Cutoff == 2')
hold on
subplot(2,2,2)
title('Cutoff == 3')
hold on
%}
if symmetric_decomposition == 1
    semilogx(sepX,P_difference1, '-rs');
    hold on
    semilogx(sepX,P_difference2, '-gs');
    hold on
    semilogx(sepX,P_difference8, '-bs');
    legend('cutoff == 1', 'cutoff == 2', 'cutoff == 8')
    legend('show')
    title('Differences between P and P_a_p_p_r_o_x with symmetric charge distribution')
    xlabel('Distance')
    ylabel('Error')
    hold off
else
    semilogx(sepX,P_difference1, '-rs');
    hold on
    semilogx(sepX,P_difference2, '-gs');
    hold on
    semilogx(sepX,P_difference8, '-bs');
    legend('cutoff == 1', 'cutoff == 2', 'cutoff == 8')
    legend('show')
    title('Differences between P and P_a_p_p_r_o_x with asymmetric charge distribution')
    xlabel('Distance')
    ylabel('Error')
    hold off
end

