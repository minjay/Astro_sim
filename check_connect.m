function flag = check_connect(k1, k2, init_sets, adj_mat)

flag = 0;
for i = init_sets{k1}
    for j = init_sets{k2}
        if adj_mat(i, j)
            flag = 1;
            return
        end
    end
end

end