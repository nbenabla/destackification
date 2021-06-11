struct X
    behavior::String
end

function makeX(first::Int64, second::Int64)
    return(X(string(first+second)))
end
