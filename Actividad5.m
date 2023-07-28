%Seminario de Inteligencia Artificial I

%Octave 5.2.0

clear all %Limpiar datos

functionNumber = input("Funcion a optimizar: "); %Numero de funcion a ejecutar

%Funcion Seleccion
function [n] = Selection(apt)
  totalAptitud = sum(apt); %Suma de las aptitudes
  aptitudSize = numel(apt); %Tama�o de las aptitudes
  
  random = rand(); %Numero aleatorio
  selectionSum = 0; %Actual de la suma
  
  for i=1:aptitudSize
    selectionSum = selectionSum + apt(i)/totalAptitud; %Sumar la aptitud actual
    
    if selectionSum >= random
      n = i; %Aptitud seleccionada
      return;
    endif
    
  end
  
  n = aptitudSize; %Seleccionar la ultima aptitud
end

%Funcion de Griewank
function [y] = griewank(x)
  d = length(x);
  sum = 0;
  prod = 1;

  for i=1:d
    xi = x(i);
    sum = sum + xi^2/4000;
    prod = prod * cos(xi/sqrt(i));
  end

  y = sum - prod + 1;

end

%Funcion de Rastrigin
function [y] = rastr(x)
  d = length(x);
  sum = 0;
  for i=1:d
    xi = x(i);
    sum = sum + (xi^2 - 10*cos(2*pi*xi));
  end

  y = 10*d + sum;

end

#Funcion de la Esfera
function [y] = spheref(x)
  d = length(x);
  sum = 0;
  for i=1:d
    xi = x(i);
    sum = sum + xi^2;
  end

  y = sum;

end

%Funcion de Dixon Price
function [y] = dixonpr(x)
  x1 = x(1);
  d = length(x);
  term1 = (x1-1)^2;

  sum = 0;
  for i=2:d
    xi = x(i);
    xold = x(i-1);
    new = i * (2*xi^2 - xold)^2;
    sum = sum + new;
  end

  y = term1 + sum;

end

G = 100; %Numero de generaciones
N = 30; %Numero de abejas
D = 2; %Tama�o de dimension

L = 35; %Numero de intentos
Pf = 10; %Numero de abejas trabajadoras
Po = N-Pf; %Numero de abejas observadoras

hold on %Multiples datos en una gr�fica
grid on %Mostrar las lineas de la gr�fica

%Seleccionar funci�n Artificial Bee Colony
switch functionNumber
  %Funci�n 1 Griewank
  case 1
    xl = -5; %Parte baja en x
    xu = 5; %Parte alta en x
    yl = -5; %Parte baja en y
    yu = 5; %Parte alta en y
    
    [x,y] = meshgrid(xl:.2:xu,yl:.2:yu); %Coordenadas a graficar

    z = []; %Coordenadas evaluadas
    iterations = (xu - xl)/.2; %Numero de coordenadas
    
    %Evaluar en funcion
    for i=1:iterations+1
      for j=1:iterations+1
        functionX = [x(i,j) y(i,j)];
        z(i,j) = griewank(functionX);
      end
    end
    
    contour(x,y,z) %Graficar
    
    xl = [-5 -5]'; %Valor bajo en x
    xu = [5 5]'; %Valor alto en x
    
    x = zeros(D,Pf); %Inicializar en cero las coordenadas
    l = zeros(1,Pf); %Inicializar en cero los intentos
    aptitud = zeros(1,Pf); %Inicializar en celo la aptitud
    fitness = zeros(1,Pf); %Inicializar en cero las evaluaciones

    %Inicializar aleatoriamente los valores
    for i=1:Pf
      x(:,i) = xl+(xu-xl).*rand(D,1);
      functionX = [x(1,i), x(2,i)];
      fitness(i) = griewank(functionX);
      
      %Calcular aptitud
      if fitness(i) >= 0
        aptitud(i) = 1/(1+fitness(i));
      else
        aptitud(i) = 1 + abs(fitness(i));
      end
    end
    
    %Generaciones
    for g=1:G
      %Abejas Empleadas
      for i=1:Pf
        %Obtener fuente de alimento
        k = randi([1, Pf]);
        
        while k == i
          k = randi([1, Pf]);
        end
        
        j = randi([1, D]); %Valor aleatorio entre la dimension del problema
        phi = 2*rand()-1; %Valor aleatorio entre -1 y 1
        
        v = x(:,i); %Respaldar vector
        v(j) = x(j,i) + phi*(x(j,i)-x(j,k));
        
        %Evaluar en funcion
        functionV = [v(1), v(2)];
        currentFunctionV = griewank(functionV);
        
        functionX = [x(1,i), x(2,i)];
        currentFunctionX = griewank(functionX);
        
        %Asignar intentos y nuevo valor
        if currentFunctionV < currentFunctionX
          x(:,i) = v;
          fitness(i) = currentFunctionV;
          l(i) = 0;
        else
          l(i) = l(i) + 1;
        end
        
      end
      
      %Calcular aptitudes
      for i=1:Pf
        if fitness(i) >= 0
          aptitud(i) = 1/(1+fitness(i));
        else
          aptitud(i) = 1 + abs(fitness(i));
        end
      end
      
      %Abejas Observadoras
      for i=1:Po
        m = Selection(aptitud); %Abeja trabajadora
        
        %Obtener fuente de alimento
        k = randi([1, Pf]);
        
        while k == m
          k = randi([1, Pf]);
        end
        
        j = randi([1, D]); %Valor aleatorio entre la dimension del problema
        phi = 2*rand()-1; %Valor aleatorio entre -1 y 1
        
        v = x(:,m); %Respaldar vector
        v(j) = x(j,m) + phi*(x(j,m)-x(j,k));
        
        %Evaluar en funcion
        functionV = [v(1), v(2)];
        currentFunctionV = griewank(functionV);
        
        functionX = [x(1,m), x(2,m)];
        currentFunctionX = griewank(functionX);
        
        %Asignar intentos y nuevo valor
        if currentFunctionV < currentFunctionX
          x(:,m) = v;
          fitness(m) = currentFunctionV;
          l(m) = 0;
        else
          l(m) = l(m) + 1;
        end
        
      end
      
      %Abejas Exploradoras
      for i=1:Pf
        %Buscar otra posicion si alcanzo el maximo de intentos
        if l(i) > L
          functionX = [x(1,i), x(2,i)];
          fitness(i) = griewank(functionX);
          x(:,i) = xl+(xu-xl).*rand(D,1);
          l(i) = 0;
        end
        
      end
      
    end
  
    [~,bestPosition] = min(fitness); %Obtener mejor posicion
    
    %Graficar minimo
    plot(x(1,bestPosition),x(2,bestPosition),'rx','LineWidth',2,'markerSize',10) %Mostrar el optimo
    
    %Mostrar resultados
    disp("Minimo");
    disp(x(1,bestPosition));
    disp(x(2,bestPosition));
    functionX = [x(1,bestPosition), x(2,bestPosition)];
    disp(griewank(functionX));
      
  %Funci�n 2 Rastrigin
  case 2
    xl = -5.12; %Parte baja en x
    xu = 5.12; %Parte alta en x
    yl = -5.12; %Parte baja en y
    yu = 5.12; %Parte alta en y
    
    [x,y] = meshgrid(xl:.2:xu,yl:.2:yu); %Coordenadas a graficar

    z = []; %Coordenadas evaluadas
    iterations = (xu - xl)/.2; %Numero de coordenadas
    
    %Evaluar en funcion
    for i=1:iterations+1
      for j=1:iterations+1
        functionX = [x(i,j) y(i,j)];
        z(i,j) = rastr(functionX);
      end
    end
    
    contour(x,y,z) %Graficar
    
    xl = [-5.12 -5.12]'; %Valor bajo en x
    xu = [5.12 5.12]'; %Valor alto en x
    
    x = zeros(D,Pf); %Inicializar en cero las coordenadas
    l = zeros(1,Pf); %Inicializar en cero los intentos
    aptitud = zeros(1,Pf); %Inicializar en celo la aptitud
    fitness = zeros(1,Pf); %Inicializar en cero las evaluaciones

    %Inicializar aleatoriamente los valores
    for i=1:Pf
      x(:,i) = xl+(xu-xl).*rand(D,1);
      functionX = [x(1,i), x(2,i)];
      fitness(i) = rastr(functionX);
      
      %Calcular aptitud
      if fitness(i) >= 0
        aptitud(i) = 1/(1+fitness(i));
      else
        aptitud(i) = 1 + abs(fitness(i));
      end
    end
    
    %Generaciones
    for g=1:G
      %Abejas Empleadas
      for i=1:Pf
        %Obtener fuente de alimento
        k = randi([1, Pf]);
        
        while k == i
          k = randi([1, Pf]);
        end
        
        j = randi([1, D]); %Valor aleatorio entre la dimension del problema
        phi = 2*rand()-1; %Valor aleatorio entre -1 y 1
        
        v = x(:,i); %Respaldar vector
        v(j) = x(j,i) + phi*(x(j,i)-x(j,k));
        
        %Evaluar en funcion
        functionV = [v(1), v(2)];
        currentFunctionV = rastr(functionV);
        
        functionX = [x(1,i), x(2,i)];
        currentFunctionX = rastr(functionX);
        
        %Asignar intentos y nuevo valor
        if currentFunctionV < currentFunctionX
          x(:,i) = v;
          fitness(i) = currentFunctionV;
          l(i) = 0;
        else
          l(i) = l(i) + 1;
        end
        
      end
      
      %Calcular aptitudes
      for i=1:Pf
        if fitness(i) >= 0
          aptitud(i) = 1/(1+fitness(i));
        else
          aptitud(i) = 1 + abs(fitness(i));
        end
      end
      
      %Abejas Observadoras
      for i=1:Po
        m = Selection(aptitud); %Abeja trabajadora
        
        %Obtener fuente de alimento
        k = randi([1, Pf]);
        
        while k == m
          k = randi([1, Pf]);
        end
        
        j = randi([1, D]); %Valor aleatorio entre la dimension del problema
        phi = 2*rand()-1; %Valor aleatorio entre -1 y 1
        
        v = x(:,m); %Respaldar vector
        v(j) = x(j,m) + phi*(x(j,m)-x(j,k));
        
        %Evaluar en funcion
        functionV = [v(1), v(2)];
        currentFunctionV = rastr(functionV);
        
        functionX = [x(1,m), x(2,m)];
        currentFunctionX = rastr(functionX);
        
        %Asignar intentos y nuevo valor
        if currentFunctionV < currentFunctionX
          x(:,m) = v;
          fitness(m) = currentFunctionV;
          l(m) = 0;
        else
          l(m) = l(m) + 1;
        end
        
      end
      
      %Abejas Exploradoras
      for i=1:Pf
        %Buscar otra posicion si alcanzo el maximo de intentos
        if l(i) > L
          functionX = [x(1,i), x(2,i)];
          fitness(i) = rastr(functionX);
          x(:,i) = xl+(xu-xl).*rand(D,1);
          l(i) = 0;
        end
        
      end
      
    end
  
    [~,bestPosition] = min(fitness); %Obtener mejor posicion
    
    %Graficar minimo
    plot(x(1,bestPosition),x(2,bestPosition),'rx','LineWidth',2,'markerSize',10) %Mostrar el optimo
    
    %Mostrar resultados
    disp("Minimo");
    disp(x(1,bestPosition));
    disp(x(2,bestPosition));
    functionX = [x(1,bestPosition), x(2,bestPosition)];
    disp(rastr(functionX));
  
  %Funcion 3 Sphere
  case 3
    xl = -5.12; %Parte baja en x
    xu = 5.12; %Parte alta en x
    yl = -5.12; %Parte baja en y
    yu = 5.12; %Parte alta en y
    
    [x,y] = meshgrid(xl:.2:xu,yl:.2:yu); %Coordenadas a graficar

    z = []; %Coordenadas evaluadas
    iterations = (xu - xl)/.2; %Numero de coordenadas
    
    %Evaluar en funcion
    for i=1:iterations+1
      for j=1:iterations+1
        functionX = [x(i,j) y(i,j)];
        z(i,j) = spheref(functionX);
      end
    end
    
    contour(x,y,z) %Graficar
    
    xl = [-5.12 -5.12]'; %Valor bajo en x
    xu = [5.12 5.12]'; %Valor alto en x
    
    x = zeros(D,Pf); %Inicializar en cero las coordenadas
    l = zeros(1,Pf); %Inicializar en cero los intentos
    aptitud = zeros(1,Pf); %Inicializar en celo la aptitud
    fitness = zeros(1,Pf); %Inicializar en cero las evaluaciones

    %Inicializar aleatoriamente los valores
    for i=1:Pf
      x(:,i) = xl+(xu-xl).*rand(D,1);
      functionX = [x(1,i), x(2,i)];
      fitness(i) = spheref(functionX);
      
      %Calcular aptitud
      if fitness(i) >= 0
        aptitud(i) = 1/(1+fitness(i));
      else
        aptitud(i) = 1 + abs(fitness(i));
      end
    end
    
    %Generaciones
    for g=1:G
      %Abejas Empleadas
      for i=1:Pf
        %Obtener fuente de alimento
        k = randi([1, Pf]);
        
        while k == i
          k = randi([1, Pf]);
        end
        
        j = randi([1, D]); %Valor aleatorio entre la dimension del problema
        phi = 2*rand()-1; %Valor aleatorio entre -1 y 1
        
        v = x(:,i); %Respaldar vector
        v(j) = x(j,i) + phi*(x(j,i)-x(j,k));
        
        %Evaluar en funcion
        functionV = [v(1), v(2)];
        currentFunctionV = spheref(functionV);
        
        functionX = [x(1,i), x(2,i)];
        currentFunctionX = spheref(functionX);
        
        %Asignar intentos y nuevo valor
        if currentFunctionV < currentFunctionX
          x(:,i) = v;
          fitness(i) = currentFunctionV;
          l(i) = 0;
        else
          l(i) = l(i) + 1;
        end
        
      end
      
      %Calcular aptitudes
      for i=1:Pf
        if fitness(i) >= 0
          aptitud(i) = 1/(1+fitness(i));
        else
          aptitud(i) = 1 + abs(fitness(i));
        end
      end
      
      %Abejas Observadoras
      for i=1:Po
        m = Selection(aptitud); %Abeja trabajadora
        
        %Obtener fuente de alimento
        k = randi([1, Pf]);
        
        while k == m
          k = randi([1, Pf]);
        end
        
        j = randi([1, D]); %Valor aleatorio entre la dimension del problema
        phi = 2*rand()-1; %Valor aleatorio entre -1 y 1
        
        v = x(:,m); %Respaldar vector
        v(j) = x(j,m) + phi*(x(j,m)-x(j,k));
        
        %Evaluar en funcion
        functionV = [v(1), v(2)];
        currentFunctionV = spheref(functionV);
        
        functionX = [x(1,m), x(2,m)];
        currentFunctionX = spheref(functionX);
        
        %Asignar intentos y nuevo valor
        if currentFunctionV < currentFunctionX
          x(:,m) = v;
          fitness(m) = currentFunctionV;
          l(m) = 0;
        else
          l(m) = l(m) + 1;
        end
        
      end
      
      %Abejas Exploradoras
      for i=1:Pf
        %Buscar otra posicion si alcanzo el maximo de intentos
        if l(i) > L
          functionX = [x(1,i), x(2,i)];
          fitness(i) = spheref(functionX);
          x(:,i) = xl+(xu-xl).*rand(D,1);
          l(i) = 0;
        end
        
      end
      
    end
  
    [~,bestPosition] = min(fitness); %Obtener mejor posicion
    
    %Graficar minimo
    plot(x(1,bestPosition),x(2,bestPosition),'rx','LineWidth',2,'markerSize',10) %Mostrar el optimo
    
    %Mostrar resultados
    disp("Minimo");
    disp(x(1,bestPosition));
    disp(x(2,bestPosition));
    functionX = [x(1,bestPosition), x(2,bestPosition)];
    disp(spheref(functionX));
    
  %Funcion 4 Dixon Price
  case 4
    xl = -10; %Parte baja en x
    xu = 10; %Parte alta en x
    yl = -10; %Parte baja en y
    yu = 10; %Parte alta en y
    
    [x,y] = meshgrid(xl:.2:xu,yl:.2:yu); %Coordenadas a graficar

    z = []; %Coordenadas evaluadas
    iterations = (xu - xl)/.2; %Numero de coordenadas
    
    %Evaluar en funcion
    for i=1:iterations+1
      for j=1:iterations+1
        functionX = [x(i,j) y(i,j)];
        z(i,j) = dixonpr(functionX);
      end
    end
    
    contour(x,y,z) %Graficar
    
    xl = [-10 -10]'; %Valor bajo en x
    xu = [10 10]'; %Valor alto en x
    
    x = zeros(D,Pf); %Inicializar en cero las coordenadas
    l = zeros(1,Pf); %Inicializar en cero los intentos
    aptitud = zeros(1,Pf); %Inicializar en celo la aptitud
    fitness = zeros(1,Pf); %Inicializar en cero las evaluaciones

    %Inicializar aleatoriamente los valores
    for i=1:Pf
      x(:,i) = xl+(xu-xl).*rand(D,1);
      functionX = [x(1,i), x(2,i)];
      fitness(i) = dixonpr(functionX);
      
      %Calcular aptitud
      if fitness(i) >= 0
        aptitud(i) = 1/(1+fitness(i));
      else
        aptitud(i) = 1 + abs(fitness(i));
      end
    end
    
    %Generaciones
    for g=1:G
      %Abejas Empleadas
      for i=1:Pf
        %Obtener fuente de alimento
        k = randi([1, Pf]);
        
        while k == i
          k = randi([1, Pf]);
        end
        
        j = randi([1, D]); %Valor aleatorio entre la dimension del problema
        phi = 2*rand()-1; %Valor aleatorio entre -1 y 1
        
        v = x(:,i); %Respaldar vector
        v(j) = x(j,i) + phi*(x(j,i)-x(j,k));
        
        %Evaluar en funcion
        functionV = [v(1), v(2)];
        currentFunctionV = dixonpr(functionV);
        
        functionX = [x(1,i), x(2,i)];
        currentFunctionX = dixonpr(functionX);
        
        %Asignar intentos y nuevo valor
        if currentFunctionV < currentFunctionX
          x(:,i) = v;
          fitness(i) = currentFunctionV;
          l(i) = 0;
        else
          l(i) = l(i) + 1;
        end
        
      end
      
      %Calcular aptitudes
      for i=1:Pf
        if fitness(i) >= 0
          aptitud(i) = 1/(1+fitness(i));
        else
          aptitud(i) = 1 + abs(fitness(i));
        end
      end
      
      %Abejas Observadoras
      for i=1:Po
        m = Selection(aptitud); %Abeja trabajadora
        
        %Obtener fuente de alimento
        k = randi([1, Pf]);
        
        while k == m
          k = randi([1, Pf]);
        end
        
        j = randi([1, D]); %Valor aleatorio entre la dimension del problema
        phi = 2*rand()-1; %Valor aleatorio entre -1 y 1
        
        v = x(:,m); %Respaldar vector
        v(j) = x(j,m) + phi*(x(j,m)-x(j,k));
        
        %Evaluar en funcion
        functionV = [v(1), v(2)];
        currentFunctionV = dixonpr(functionV);
        
        functionX = [x(1,m), x(2,m)];
        currentFunctionX = dixonpr(functionX);
        
        %Asignar intentos y nuevo valor
        if currentFunctionV < currentFunctionX
          x(:,m) = v;
          fitness(m) = currentFunctionV;
          l(m) = 0;
        else
          l(m) = l(m) + 1;
        end
        
      end
      
      %Abejas Exploradoras
      for i=1:Pf
        %Buscar otra posicion si alcanzo el maximo de intentos
        if l(i) > L
          functionX = [x(1,i), x(2,i)];
          fitness(i) = dixonpr(functionX);
          x(:,i) = xl+(xu-xl).*rand(D,1);
          l(i) = 0;
        end
        
      end
      
    end
  
    [~,bestPosition] = min(fitness); %Obtener mejor posicion
    
    %Graficar minimo
    plot(x(1,bestPosition),x(2,bestPosition),'rx','LineWidth',2,'markerSize',10) %Mostrar el optimo
    
    %Mostrar resultados
    disp("Minimo");
    disp(x(1,bestPosition));
    disp(x(2,bestPosition));
    functionX = [x(1,bestPosition), x(2,bestPosition)];
    disp(dixonpr(functionX));
    
  otherwise
    disp("Numero no valido")

endswitch