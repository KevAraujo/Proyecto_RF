clear all; close all; clc;

%% Preguntar si el usuario quiere leer o escribir un archivo Touchstone
disp('¿Qué desea hacer?');
disp('1: Leer archivo Touchstone');
disp('2: embedding archivos sNp');
disp('3: Escribir archivo Touchstone');
opcion_archivo = input('Ingrese el número de la opción: ');

if opcion_archivo == 1
 % Leer los datos del archivo
    [param,Z0] = Read_print;
    freq = [];
    parametros_S = [];
    
    freq=param.Frequencies;
    parametros_S=param.Parameters;
    N=2;

elseif opcion_archivo == 2

    % Leer los datos del archivo
    [param,Z0] = Read_print;
    freq = [];
    parametros_S = [];

    [parametros_S] = embeding(param);
    
    freq=param.Frequencies;
    %parametros_S=param.Parameters;
    N=2;



elseif opcion_archivo == 3
  Z = []; % Initialize Z to an empty matrix
[nV, nz, frec1, frec2, paso] = Solicitudes();

% Precalculate the frequency range and number of steps
frequencies = frec1:paso:frec2;
num_frequencies = length(frequencies);

% Initialize matrices
V_nodos = zeros(nV, 2);
mat_imp = zeros(nz, 1, num_frequencies);
n_matriz = 0;
matriz_parametros = zeros(nV, nV);
matriz_nodos = zeros(nz, 2);

% Initialize impedance structure
Impedancias = struct();

for i = 1:nz
while true
    soli = input('¿Qué desea ingresar? R=1; L=2; C=3: ');
    if ismember(soli, [1, 2, 3])
        break; % Exit loop if input is valid
    else
        disp('Error: No seleccionó la entrada adecuada. Inténtelo de nuevo.');
    end
end
    
%% Input component values and check for validity
if soli == 2  % Inductor
    L_value = input(['L' num2str(i) ' = ']);
    assert(L_value > 0, 'El valor de la inductancia debe ser positivo.');
    x(i) = 1i * 2 * pi * L_value;  % Inductancia en forma de reactancia
elseif soli == 3  % Capacitor
    C_value = input(['C' num2str(i) ' = ']);
    assert(C_value > 0, 'El valor de la capacitancia debe ser positivo.');
    x(i) = -1i / (2 * pi * C_value);  % Capacitancia en forma de reactancia
else  soli == 1 % Resistor
    R_value = input(['R' num2str(i) ' = ']);
    assert(R_value >= 0, 'El valor de la resistencia no puede ser negativo.');
    x(i) = R_value;  % Resistencia
end
    
    % Input nodes
    matriz_nodos(i, :) = str2double(split(input('Nodos: ', 's')));
    
    %% Calcula impedancias para cada frecuencia
    for idx = 1:num_frequencies
        frec = frequencies(idx);
        if soli == 2
            z(i) = frec * x(i);
        elseif soli == 3
            z(i) = x(i) / frec;
        else
         z(i) = x(i);
        end
        mat_imp(i, :, idx) = z(i);
        fieldName = sprintf('Z%d', i);
        Impedancias.(fieldName) = z(i);
    end
end

n_nodos = max(matriz_nodos(:)) + 1;

% Initialize numeric matrices for Y calculations
V_nodos_num = zeros(nV, 2);
Vn = zeros(n_nodos - 1, 1);

% Initialize numeric current variables
z_matriz = cell(nz, 1);
matriz_Y = zeros(n_nodos - 1, n_nodos - 1, num_frequencies);
matriz_y = matriz_Y;

% Input node voltages
for i = 1:nV
    V_nodos_num(i, :) = str2double(split(input(['Nodos de V' num2str(i) ': '], 's')));
    np = max(V_nodos_num(i, :));
    Vn(np) = 1;
end

pos_V = find(Vn == 0);

% Calculate Y matrix for each frequency
for idx = 1:num_frequencies
    frec = frequencies(idx);
    for i = 1:nz
        z_matriz{i} = zeros(n_nodos - 1);
        nmax = max(matriz_nodos(i, :));
        nmin = min(matriz_nodos(i, :));
        
        if nmin > 0
            z_matriz{i}(nmin, nmin) = 1 / mat_imp(i, :, idx);
            z_matriz{i}(nmin, nmax) = -1 / mat_imp(i, :, idx);
            z_matriz{i}(nmax, nmin) = -1 / mat_imp(i, :, idx);
            z_matriz{i}(nmax, nmax) = 1 / mat_imp(i, :, idx);
        else
            z_matriz{i}(nmax, nmax) = 1 / mat_imp(i, :, idx);
        end
        
        matriz_Y(:, :, idx) = matriz_Y(:, :, idx) + z_matriz{i};
    end
    
    % Calculate Y matrix with row/column suppression for floating nodes
    for i = 1:(length(matriz_Y(:, :, idx)) - nV)
        q = pos_V(i);
        for j = 1:length(matriz_Y(:, :, idx))
            for k = 1:length(matriz_Y(:, :, idx))
                matriz_y(j, k, idx) = matriz_Y(j, k, idx) - ...
                    (matriz_Y(j, q, idx) * matriz_Y(q, k, idx)) / matriz_Y(q, q, idx);
            end
        end
        matriz_Y(:, :, idx) = matriz_y(:, :, idx);
    end
end

% Remove suppressed nodes
if ~isempty(pos_V)
    matriz_Y(pos_V, :, :) = [];
    matriz_Y(:, pos_V, :) = [];
end

% Create ABCD matrix for each frequency
Matriz_ABCD = zeros(nV, nV, num_frequencies);
%% selección de parametros 
while true
    disp('Seleccione la matriz que desea ver:');
    disp('1: Matriz Y (admitancia)');
    disp('2: Matriz Z (impedancia)');
    disp('3: Matriz ABCD');
    disp('4: Matriz S');
    disp('0: Salir');
    
    opcion = input('Ingrese el número de la opción: ');
    
    switch opcion
        case 1
           for idx = 1:num_frequencies
            if isequal(V_nodos_num(1, :), V_nodos_num(2, :))
                disp(['Frecuencia ', num2str(frequencies(idx)), ': Caso Especial, no hay matriz Y']);
            else
                disp(['Frecuencia ', num2str(frequencies(idx)), ': Matriz de parámetros Y']);
                disp(matriz_Y(:,:,idx)); % Mostrar la matriz Y para esta frecuencia específica
            end
        end
            
        case 2
            for idx = 1:num_frequencies
                detY = det(matriz_Y(:,:,idx));
                if detY == 0
                    disp(['Frecuencia ', num2str(frequencies(idx)), ': La matriz Y no es invertible.']);
                else
                    Z(:,:,idx) = inv(matriz_Y(:,:,idx));
                    disp(['Matriz Z en frecuencia ', num2str(frequencies(idx)), ':']);
                    disp(Z(:,:,idx));
                end
            end
            
        case 3
            for idx = 1:num_frequencies
                if det(matriz_Y(:,:,idx)) == 0 && ~isequal(V_nodos_num(1, :), V_nodos_num(2, :))
                    A = 1; 
                    B = z(1); 
                    C = 0; 
                    D = 1;
                    disp(['Frecuencia ', num2str(frequencies(idx)), ': Caso especial - matriz Y no es invertible.']);
                elseif (V_nodos_num(1, :) == V_nodos_num(2, :))
            % Caso Especial 2: Los nodos de V1 y V2 son iguales
            A = 1;
            B = 0;  % Asignamos un valor específico para B cuando los nodos son iguales
            C = (z(1))^-1; % Utilizamos la impedancia calculada
            D = 1;
                else
                    Z_curr = Z(:,:,idx);
                    A = Z_curr(1, 1) / Z_curr(2, 1);
                    B = det(Z_curr) / Z_curr(2, 1);
                    C = 1 / Z_curr(2, 1);
                    D = Z_curr(2, 2) / Z_curr(2, 1);
                end
                
                Matriz_ABCD(:,:,idx) = [A, B; C, D];
                disp(['Matriz ABCD en frecuencia ', num2str(frequencies(idx)), ':']);
                disp(Matriz_ABCD(:,:,idx));
            end
       case 4
    % Solicitar impedancia de referencia Z0
    Z0 = input('Ingrese la impedancia de referencia Z0 (50 Ω, solo presione enter para valor predeterminado): ');
    if isempty(Z0)
        Z0 = 50; % Valor por defecto
    end
% Llamar a la función para convertir Z a S
Matriz_S = Z_to_S(matriz_Y,Z, Z0);

% Mostrar los resultados
for idx = 1:size(Matriz_S, 3)
    disp(['Frecuencia ', num2str(idx), ': Matriz S']);
    disp(Matriz_S(:, :, idx));
end
%% escribir archivo sNp
% Preguntar si se desea guardar en un archivo sNp
    guardar = input('¿Desea guardar los parámetros S en un archivo .sNp? (s/n): ', 's');
    if lower(guardar) == 's'
        Imprimir(Matriz_S,nV,frec1,frec2,paso);
    end

        case 0
            disp('Saliendo del programa.');
            break;
            
        otherwise
            disp('Opción no válida. Por favor seleccione 1, 2, 3 o 0 para salir.');
    end
    end

else
    disp('Opción no válida. Saliendo.');
end
%% graficar
graficarSeleccion(parametros_S, N, freq);

if opcion_archivo == 2
    guardar = input('¿Desea guardar los parámetros S en un archivo .sNp? (s/n): ', 's');
    if lower(guardar) == 's'
    frec1=freq(1);
    frec2=freq(length(freq));
    paso=freq(2)-freq(1);
    nV=2;
    Imprimir(parametros_S,nV,frec1,frec2,paso);
    end
end
