
function vecOuterProd!(C, a, b, α, β)
    @turbo for i ∈ eachindex(b)
        for j ∈ eachindex(a)
            C[j,i] = β*C[j,i] + α*a[j]*b[i]
        end
    end
end

function vecOuterProd!(C, a, b)
    @turbo for i ∈ eachindex(b)
        for j ∈ eachindex(a)
            C[j,i] = a[j]*b[i]
        end
    end
end

function vecOuterProdPlusI!(C, a, b)
    @turbo for i ∈ eachindex(b)
        for j ∈ eachindex(a)
            C[j,i] = a[j]*b[i] + I[j,i]
        end
    end
end
