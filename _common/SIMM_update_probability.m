function miu = SIMM_update_probability(k,model,estS,c,miu)
epsilon = 1e-20;%原来为1e-6，但是感觉太大，万一只有一个模型的Lk为0，则赋值后占主导

Lk_last = [];
if isempty(estS.Lk{k,1})|isempty(estS.Lk{k,2})
        miu= miu;         %若似然为空，则不对miu更新，即miu保持不变
else

    Lk = [estS.Lk{k,1}  estS.Lk{k,2} ];
    if all(Lk(:) == 0)     %全零矩阵
        Lk_min = epsilon;  % 指定默认值（如 eps）
    else
        Lk_min = min(Lk(Lk ~= 0));  %找出Lk矩阵中非零元素的最小值
    end

    for mll =1:model.M
        LLLk_d = estS.Lk{k,mll};
        if LLLk_d == 0
            LLLk_d = Lk_min*epsilon;
        end
        Lk_last = [Lk_last LLLk_d];  %存
    end

    c_norm = max(Lk_last .* c');  %归一化常数

    for mll =1:model.M
        miu(mll,:) = Lk_last(:,mll) * c(mll,:)./c_norm;
    end
    miu=miu/max(miu);

end
end