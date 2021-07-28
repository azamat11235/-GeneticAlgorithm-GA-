/// Генетический алгоритм
/// 09.03.2020
/// Султонов Азамат
///   Критерий останова: достижение априорно заданного значения целевой функции с заданной точностью
///   Отбор особей для скрещивания: турнирный
///   Кроссинговер: двухточечный
///   Мутация: изменение случайно выбранного бита

const N = 10; {размер начальной популяции}
	e = 1/16384; {точность решения} // 0,00006103515
	resultX = 3.81046;
	resultY = 30.119;
type list = ^node;
	node = record x,funcValue: real; next: list end;
var population,offspring: list; numOfIter,popSize,mode, countOfCros, countOfMut, stdOut: integer; probabilityOfMutation, probabilityOfCros: real;
	resFile: text;

function f(x: real): real; {данная функция}
	begin
		f := (x-1)*(x-1)*(x-1)*(x-1)*(x-1) * (x-0.05) * (x-3) * (x-3.5) * (1-exp(x-3.95)) * ln(x+0.22); { (x-1)^5 * (x-0.05) * (x-3) * (x-3.5) * (1-e^(x-3.95)) * ln(x+0.22) }
	end;

procedure delList(L: list);
	begin
		if L <> nil then begin
			delList(L^.next);
			dispose(L);
		end;
	end;

procedure addToList(var L, elem: list);
	var p,q: list;
	begin
		if L <> nil then begin
			if elem^.funcValue < L^.funcValue then begin
				p := L;
				while (p^.next <> nil) and (elem^.funcValue < p^.next^.funcValue) do p := p^.next;
				if p^.next = nil then begin
					new(p^.next);
					p := p^.next;
					p^.x := elem^.x;
					p^.funcValue := elem^.funcValue;
					p^.next := nil;
				end else begin
					if elem^.x <> p^.next^.x then begin
						new(q);
						q^.x := elem^.x;
						q^.funcValue := elem^.funcValue;
						q^.next := p^.next;
						p^.next := q;
					end;
				end;
			end else begin
				if elem^.x <> L^.x then begin
					new(p);
					p^.x := elem^.x;
					p^.funcValue := elem^.funcValue;
					p^.next := L;
					L := p;
				end;
			end;
		end else begin
			new(L);
			L^.x := elem^.x;
			L^.funcValue := elem^.funcValue;
			L^.next := nil;
		end;
	end;

procedure writeBin(r: real; i: integer);
	begin
		if r < 0 then begin
			write('-');
			r := abs(r);
		end;
		if r <> 0 then begin
			while r - exp((i-1)*ln(2)) < 0 do begin
				write('0');
				i := i - 1;
			end;
			if trunc(r) div 2 <> 0 then writeBin(trunc(r) div 2, i-1);
			if trunc(r) mod 2 = 1 then write('1')
			else write('0');
			r := r - trunc(r);
			if r > 0 then begin
				write(',');
				r := r * 2;
				if r >= 1 then begin
					write('1');
					r := r - 1;
				end else write('0');
				while r <> 0 do begin
					r := r * 2;
					if r >= 1 then begin
						write('1');
						r := r - 1;
					end else write('0');
				end;
			end;
		end else begin
			while i <> 0 do begin
				write('0');
				i := i - 1;
			end;
		end;
	end;

procedure writePop(var population: list; numOfIter: integer);
	var p: list; i: integer;
	begin
		p := population;
		i := 1;
		writeln('               ---------------------------------------------------------------------');
		writeln('               --    --                           Особь                           --');
		writeln('               -- №  ---------------------------------------------------------------');
		writeln('               --    -- Десятичное число --  Двоичное число  -- Приспособленность --');
		writeln('               ---------------------------------------------------------------------');
		while p <> nil do begin
			write('               -- ', i:2, ' -- ', p^.x*exp(-14*ln(2)):0:14, ' -- ');
			writeBin(p^.x, 16);
			writeln(' --', p^.funcValue:18:4, ' --');
			i := i + 1;
			p := p^.next;
		end;
		writeln('               ---------------------------------------------------------------------');
	end;

function stopCriteria(var L: list): boolean;
	var p: list;
	begin
		p := L;
		while (p <> nil) and (abs(p^.funcValue - resultY) > e) do
			p := p^.next;
		if p = nil then stopCriteria := false
		else begin
			writeln('> Результат:');
			writeln('     x = ', p^.x*exp(-14*ln(2)):0:14);
			writeln('     f(', p^.x*exp(-14*ln(2)):0:14, ') = ', p^.funcValue:0:4);
			writeln('     Число поколений: ', numOfIter);
			stopCriteria := true;
		end;
	end;

procedure createPop(var population: list); {формрование начальной популяции}
	var p: list; x: real; i: integer;
	begin
		new(population);
		p := population;
		p^.x := 0;
		p^.funcValue := f(0);
		x := 0.4;
		for i:=1 to N-1 do begin
			new(p^.next);
			p := p^.next;
			p^.x := trunc(x*exp(14*ln(2)));
			p^.funcValue := f(p^.x * exp(-14*ln(2)));
			x := x + 0.4;
		end;
		p^.next := nil;
	end;

procedure selection(var population,individ: list; var popSize: integer); {турнирный отбор}
	var p,q: list; i,j: integer;
	begin
		i := random(popSize) + 1;
		p := population;
		for j:=1 to i-1 do p := p^.next;
		j := random(popSize) + 1;
		while j = i do j := random(popSize) + 1;
		q := population;
		for i:=1 to j-1 do q := q^.next;
		if p^.funcValue > q^.funcValue then individ := p
		else individ := q;
	end;

procedure crossing(var population, offspring: list; var popSize, countOfCros: integer; probabilityOfCros: real); {двухточечный кроссинговер}
	var i,m,n: integer; individ_1,individ_2,r: list;
	begin
		if mode = 2 then begin
			writeln('                _____________');
			writeln('               | Скрещивания |');
			writeln('               |_____________|');
			writeln('               ---------------------------------------------------------------------');
			writeln('               -- No --      Особь 1     --      Особь 2     --      Потомки      --');
			writeln('               ---------------------------------------------------------------------');
		end;
		r := nil;
		for i:=1 to countOfCros do begin
			selection(population, individ_1, popSize);
			selection(population, individ_2, popSize);
			if mode = 2 then begin
				write('               -- ', i:2, ' -- ');
				writeBin(individ_1^.x, 16);
				write(' -- ');
				writeBin(individ_2^.x, 16);
				write(' -- ');
			end;
			if random(101) <= probabilityOfCros then begin
				if r = nil then begin
					new(offspring);
					r := offspring;
				end else begin
					new(r^.next);
					r := r^.next;
				end;
				n := random(14) + 1;
				m := n;
				while m = n do m := random(14) + 1;
				if m > n then begin
					m := m + n;
					n := m - n;
					m := m - n;
				end;
				r^.x := individ_1^.x + exp(n*ln(2)) * trunc(individ_1^.x * exp(-n*ln(2))) - exp(m*ln(2)) * trunc(individ_1^.x * exp(-m*ln(2))) + exp(m*ln(2)) * trunc(individ_2^.x * exp(-m*ln(2))) - exp(n*ln(2)) * trunc(individ_2^.x * exp(-n*ln(2)));
				r^.funcValue := f(r^.x * exp(-14*ln(2)));
				if mode = 2 then writeBin(r^.x, 16);
				new(r^.next);
				r := r^.next;
				r^.x := individ_2^.x + exp(n*ln(2)) * trunc(individ_2^.x * exp(-n*ln(2))) - exp(m*ln(2)) * trunc(individ_2^.x * exp(-m*ln(2))) + exp(m*ln(2)) * trunc(individ_1^.x * exp(-m*ln(2))) - exp(n*ln(2)) * trunc(individ_1^.x * exp(-n*ln(2)));
				r^.funcValue := f(r^.x * exp(-14*ln(2)));
				if mode = 2 then begin
					writeln('  --');
					write('               --    --                  --                  -- ');
					writeBin(r^.x, 16);
					writeln('  -- ');
					writeln('               ---------------------------------------------------------------------');
				end;
			end else if mode = 2 then begin
						writeln('----------------- --');
						writeln('               ---------------------------------------------------------------------');
					end;
		end;
		if mode = 2 then writeln;
		if r <> nil then r^.next := nil
		else offspring := nil;
	end;

procedure mutations(var offspring: list; countOfMut: integer; var probabilityOfMutation: real); {изменение случайно выбранного бита}
	var p: list; i,j,r,popSize: integer;
	begin
		if mode=2 then begin
			writeln('                _________');
			writeln('               | Мутации |');
			writeln('               |_________|');
		end;
		p := offspring;
		popSize := 0;
		while p <> nil do begin
			popSize := popSize + 1;
			p := p^.next;
		end;
		if popSize <> 0 then begin
			if mode = 2 then begin
				writeln('               ---------------------------------------------------------------------');
				writeln('               -- No --      Особь       -- Результат мутации -- Мутировавший ген --');
				writeln('               ---------------------------------------------------------------------');
			end;
			for i:=1 to countOfMut do begin
				r := random(popSize) + 1;
				p := offspring;
				for j:=1 to r-1 do p := p^.next;
				if mode = 2 then begin
					write('               -- ', i:2, ' -- ');
					writeBin(p^.x, 16);
					write(' -- ');
				end;
				if random(101) <= probabilityOfMutation then begin
					r := random(15);
					j := trunc(exp(r*ln(2)));
					if j mod 2 <> 0 then j := j + 1;
					if (trunc(p^.x) and j) = 0 then
						p^.x := p^.x + j
					else p^.x := p^.x - j;
					p^.funcValue := f(p^.x * exp(-14*ln(2)));
					if mode = 2 then begin
						writeBin(p^.x, 16);
						writeln('  --', r+1:10, '        -- ');
					end;
				end else if mode = 2 then writeln('----------------- -- ---------------- --');
			end;
		end;
		if mode = 2 then writeln('               ---------------------------------------------------------------------');
	end;

procedure selection2(var population, offspring: list); {отбор наиболее приспособленных среди родителей и потомков}
	var newPop,p: list; i: integer;
	begin
		newPop := nil;
		p := population;
		while p <> nil do begin
			addToList(newPop,p);
			p := p^.next;
		end;
		p := offspring;
		while p <> nil do begin
			addToList(newPop,p);
			p := p^.next;
		end;
		p := newPop;
		for i:=1 to N-1 do p := p^.next;
		if p <> nil then begin
			delList(p^.next);
			p^.next := nil;
		end;
		delList(population);
		delList(offspring);
		population := newPop;
	end;

begin
	randomize;
	popSize := N;
	createPop(population);
	write('> Режим программы (1 - основной, 2 - тестовый): ');
	readln(mode);
	write('> Вывод данных (1 - в терминал, 2 - в файл): ');
	readln(stdOut);
	write('> Количество скрещиваний: ');
	readln(countOfCros);
	write('> Вероятность скрещивания (%): ');
	readln(probabilityOfCros);
	write('> Количество мутаций: ');
	readln(countOfMut);
	write('> Вероятность мутации (%): ');
	readln(probabilityOfMutation);
	numOfIter := 1;
	if stdOut = 2 then begin
		assign(resFile, 'resFile.txt');
		rewrite(resFile);
		OUTPUT := resFile;
	end;
	if mode = 2 then begin
		writeln('                                    ===============================');
		writeln('                                             ПОКОЛЕНИЕ № ', numOfIter);
		writeln('                                    ===============================');
		writeln;
		writePop(population, numOfIter);
	end;
	while not stopCriteria(population) do begin
		numOfIter := numOfIter + 1;
		if mode = 2 then begin
			writeln; writeln;
			writeln('                                    ===============================');
			writeln('                                             ПОКОЛЕНИЕ № ', numOfIter);
			writeln('                                    ===============================');
		end;
		crossing(population, offspring, popSize, countOfCros, probabilityOfCros);
		mutations(offspring, countOfMut, probabilityOfMutation);
		if mode = 2 then writeln;
		selection2(population, offspring);
		if mode = 2 then begin
			writeln('                _________________');
			writeln('               | Новая популяция |');
			writeln('               |_________________|');
			writePop(population, numOfIter);
		end;
	end;
	//if stdOut = 2 then close(resFile);
	delList(population);
end.
