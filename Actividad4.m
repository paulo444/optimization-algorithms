%Seminario de Inteligencia Artificial I

%Octave 5.2.0

clear all %Limpiar datos

functionNumber = input("Funcion a optimizar: "); %Numero de funcion a ejecutar

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

G = 150; %Numero de generaciones
N = 50; %Numero de elementos
D = 2; %Tama�o de dimension

F = 1.2; %Factor de aplicaci�n
CR = 0.6; %Constante de recombinaci�n

hold on %Multiples datos en una gr�fica
grid on %Mostrar las lineas de la gr�fica

%Seleccionar funci�n Evoluci�n Diferencial
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
    
    x = zeros(D,N); %Inicializar en cero las coordenadas
    fitness = zeros(1,N); %Inicializar en cero las evaluaciones

    %Inicializar aleatoriamente los valores
    for i=1:N
      x(:,i) = xl+(xu-xl).*rand(D,1);
      functionX = [x(1,i), x(2,i)];
      fitness(i) = griewank(functionX);
    end
    
    %Generaciones
    for g=1:G
      %Evaluar nuevas posiciones
      for i=1:N
        %Obtener tres posiciones diferentes entre si y de i
        r1 = randi([1, N]);
        
        while r1 == i
          r1 = randi([1, N]);
        end
        
        r2 = randi([1, N]);
        
        while r2 == r1 || r2 == i
          r2 = randi([1, N]);
        end
        
        r3 = randi([1, N]);
        
        while r3 == r2 || r3 == r1 || r3 == i
          r3 = randi([1, N]);
        end
        
        v = x(:,r1) + F*(x(:,r2)-x(:,r3)); %Vector mutado
        
        u = zeros(D,1); %Vector de prueba
        
        %Recombinacion
        for j=1:D
          ra = rand();
          if ra <= CR
            u(j) = v(j);
          else
            u(j) = x(j,i);
          end
        end
        
        %Seleccion
        functionTest = [u(1), u(2)];
        fTest = griewank(functionTest);
        
        if fTest < fitness(i)
          x(:,i) = [u(1), u(2)]';
          fitness(i) = fTest;
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
    
    x = zeros(D,N); %Inicializar en cero las coordenadas
    fitness = zeros(1,N); %Inicializar en cero las evaluaciones

    %Inicializar aleatoriamente los valores
    for i=1:N
      x(:,i) = xl+(xu-xl).*rand(D,1);
      functionX = [x(1,i), x(2,i)];
      fitness(i) = rastr(functionX);
    end
    
    %Generaciones
    for g=1:G
      %Evaluar nuevas posiciones
      for i=1:N
        %Obtener tres posiciones diferentes entre si y de i
        r1 = randi([1, N]);
        
        while r1 == i
          r1 = randi([1, N]);
        end
        
        r2 = randi([1, N]);
        
        while r2 == r1 || r2 == i
          r2 = randi([1, N]);
        end
        
        r3 = randi([1, N]);
        
        while r3 == r2 || r3 == r1 || r3 == i
          r3 = randi([1, N]);
        end
        
        v = x(:,r1) + F*(x(:,r2)-x(:,r3)); %Vector mutado
        
        u = zeros(D,1); %Vector de prueba
        
        %Recombinacion
        for j=1:D
          ra = rand();
          if ra <= CR
            u(j) = v(j);
          else
            u(j) = x(j,i);
          end
        end
        
        %Seleccion
        functionTest = [u(1), u(2)];
        fTest = rastr(functionTest);
        
        if fTest < fitness(i)
          x(:,i) = [u(1), u(2)]';
          fitness(i) = fTest;
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
    
    x = zeros(D,N); %Inicializar en cero las coordenadas
    fitness = zeros(1,N); %Inicializar en cero las evaluaciones

    %Inicializar aleatoriamente los valores
    for i=1:N
      x(:,i) = xl+(xu-xl).*rand(D,1);
      functionX = [x(1,i), x(2,i)];
      fitness(i) = spheref(functionX);
    end
    
    %Generaciones
    for g=1:G
      %Evaluar nuevas posiciones
      for i=1:N
        %Obtener tres posiciones diferentes entre si y de i
        r1 = randi([1, N]);
        
        while r1 == i
          r1 = randi([1, N]);
        end
        
        r2 = randi([1, N]);
        
        while r2 == r1 || r2 == i
          r2 = randi([1, N]);
        end
        
        r3 = randi([1, N]);
        
        while r3 == r2 || r3 == r1 || r3 == i
          r3 = randi([1, N]);
        end
        
        v = x(:,r1) + F*(x(:,r2)-x(:,r3)); %Vector mutado
        
        u = zeros(D,1); %Vector de prueba
        
        %Recombinacion
        for j=1:D
          ra = rand();
          if ra <= CR
            u(j) = v(j);
          else
            u(j) = x(j,i);
          end
        end
        
        %Seleccion
        functionTest = [u(1), u(2)];
        fTest = spheref(functionTest);
        
        if fTest < fitness(i)
          x(:,i) = [u(1), u(2)]';
          fitness(i) = fTest;
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
    
    x = zeros(D,N); %Inicializar en cero las coordenadas
    fitness = zeros(1,N); %Inicializar en cero las evaluaciones

    %Inicializar aleatoriamente los valores
    for i=1:N
      x(:,i) = xl+(xu-xl).*rand(D,1);
      functionX = [x(1,i), x(2,i)];
      fitness(i) = dixonpr(functionX);
    end
    
    %Generaciones
    for g=1:G
      %Evaluar nuevas posiciones
      for i=1:N
        %Obtener tres posiciones diferentes entre si y de i
        r1 = randi([1, N]);
        
        while r1 == i
          r1 = randi([1, N]);
        end
        
        r2 = randi([1, N]);
        
        while r2 == r1 || r2 == i
          r2 = randi([1, N]);
        end
        
        r3 = randi([1, N]);
        
        while r3 == r2 || r3 == r1 || r3 == i
          r3 = randi([1, N]);
        end
        
        v = x(:,r1) + F*(x(:,r2)-x(:,r3)); %Vector mutado
        
        u = zeros(D,1); %Vector de prueba
        
        %Recombinacion
        for j=1:D
          ra = rand();
          if ra <= CR
            u(j) = v(j);
          else
            u(j) = x(j,i);
          end
        end
        
        %Seleccion
        functionTest = [u(1), u(2)];
        fTest = dixonpr(functionTest);
        
        if fTest < fitness(i)
          x(:,i) = [u(1), u(2)]';
          fitness(i) = fTest;
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