function Matriz_S = Z_to_S(matriz_Y,Z, Z0)
% Z_to_S - Convierte parámetros Z (impedancia) a parámetros S (dispersión)
% 
% Sintaxis:
%   Matriz_S = Z_to_S(Z, Z0)
%
% Entradas:
%   Z  - Matriz de impedancia de tamaño (n x n x f), donde:
%        n: número de puertos.
%        f: número de frecuencias.
%   Z0 - Impedancia de referencia (por defecto 50 Ω si no se proporciona).
%
% Salidas:
%   Matriz_S - Matriz de parámetros S del mismo tamaño que Z.

    % Validar la entrada Z0 (impedancia de referencia)
    if nargin < 2 || isempty(Z0)
        Z0 = 50; % Valor por defecto
    end
    
    % Validar dimensiones de la matriz Z
    [n, m, f] = size(Z);
    if n ~= m
        error('La matriz Z debe ser cuadrada en cada frecuencia.');
    end
    
    % Inicializar matriz S
    Matriz_S = zeros(n, n, f);
    I = eye(n); % Matriz identidad
    
    % Convertir Z a S para cada frecuencia
    for idx = 1:f
        Z_curr = Z(:, :, idx);
        
        % Evitar errores por singularidad
        S_den = Z_curr + Z0 * I;
        if rcond(S_den) < 1e-12
            warning(['Matriz Z+Z0I es singular o cercana a serlo en frecuencia índice ', num2str(idx)]);
            continue;
        end
        
        % Fórmula de conversión
        S_num = Z_curr - Z0 * I; % Numerador
        Matriz_S(:, :, idx) = S_num / S_den; % División matricial
    end
    if isempty(Matriz_S)
        [n, m, f] = size(matriz_Y);
        n_matriz=0; % Reiniciar cuenta
        for frec=1:1:f
        n_matriz=n_matriz+1;

        b=(1/Z0)*speye(length(matriz_Y(:,:,n_matriz)))-matriz_Y(:,:,n_matriz);
        a=matriz_Y(:,:,n_matriz)+(1/Z0)*speye(length(matriz_Y(:,:,n_matriz)));
        Matriz_S(:,:,n_matriz)=b/a;
        end
    end
end
