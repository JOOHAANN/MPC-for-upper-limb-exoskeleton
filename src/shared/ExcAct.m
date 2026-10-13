function Activation = ExcAct(Excitation, A)
    Activation = (exp(A * Excitation) - 1) / (exp(A) - 1);
end
