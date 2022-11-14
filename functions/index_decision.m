function id = index_decision(nu,n,h)

    sum_nu = 0;
    for i = 1:n
        id_start(i,1) = sum_nu + 1;
        id_end(i,1) = sum_nu + nu(i)*h;
        sum_nu = sum_nu + nu(i)*h;
    end
    id =[id_start, id_end];
end