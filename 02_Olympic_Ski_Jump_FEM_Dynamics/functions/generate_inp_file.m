function [elements_output] = generate_inp_file(filename, fundamental_nodes, connections, isolated_nodes, properties_list, target_prop_id)
    % --- INTESTAZIONE E DEFINIZIONE INPUT ---
    % fundamental_nodes: Matrice dei nodi principali [ID_nodo, vincolo_X, vincolo_Y, vincolo_Rot, coord_X, coord_Y]
    % connections: Cell array che definisce come i nodi sono collegati {ID_nodo_iniziale, ID_nodo_finale, 'line'/'arc', dati_extra, Lunghezza_target_elemento, ID_proprietà}
    % isolated_nodes: Matrice di nodi isolati (stesso formato dei fundamental_nodes) che non fanno parte di linee continue
    % properties_list: Matrice delle proprietà dei materiali/sezioni [ID_prop, massa_lineica, rigidezza_assiale_EA, rigidezza_flessionale_EJ]
    % target_prop_id: L'ID della proprietà di cui si vogliono estrarre gli elementi finali discretizzati.
    %
    % OUTPUT:
    % elements_output: Vettore colonna contenente la numerazione finale degli elementi trave che possiedono la proprietà target_prop_id.
    
    % Apre (o crea) il file di testo in modalità scrittura ('w'). 'fid' è l'identificativo del file.
    fid = fopen(filename, 'w');
    if fid == -1
        % Se fid è -1 significa che c'è stato un errore (es. permessi negati o disco pieno), quindi interrompe lo script.
        error('Impossibile aprire o creare il file %s', filename);
    end
    
    % ==========================================
    % --- 1. MAPPATURA NODI (RE-INDEXING) ---
    % ==========================================
    % Unisce i nodi fondamentali e quelli isolati in un'unica matrice.
    original_points = [fundamental_nodes; isolated_nodes];
    
    % Calcola quanti nodi principali "fissi" ci sono in totale.
    num_fixed = size(original_points, 1);
    
    % Crea un vettore 'lookup' che farà da "dizionario" per tradurre i vecchi ID nei nuovi ID continui.
    % La dimensione è pari all'ID massimo presente tra i nodi originali.
    lookup = zeros(max(original_points(:,1)), 1);
    
    % Inizializza la matrice che conterrà le coordinate e i vincoli dei nodi finali riordinati.
    final_nodes = zeros(num_fixed, 6);
    
    for i = 1:num_fixed
        % Estrae l'ID originale del nodo.
        old_id = original_points(i, 1);
        % Il nuovo ID è semplicemente l'indice del ciclo 'i' (da 1 a num_fixed).
        new_id = i; 
        % Salva nel "dizionario" la corrispondenza: se cerco old_id, mi restituisce new_id.
        lookup(old_id) = new_id;
        % Popola la riga del nuovo nodo mantenendo i vincoli e le coordinate (colonne da 2 a 6) originali.
        final_nodes(new_id, :) = [new_id, original_points(i, 2:6)];
    end
    
    % Imposta il contatore per i futuri nodi interpolati (che verranno creati durante la discretizzazione).
    next_id = num_fixed + 1; 
    
    % Inizializza la matrice che conterrà la topologia degli elementi finiti [ID_nodo1, ID_nodo2, ID_proprietà].
    final_beams = []; 
    
    % ==========================================================
    % --- 2. DISCRETIZZAZIONE E ASSEGNAZIONE PROPRIETÀ ---
    % ==========================================================
    for i = 1:size(connections, 1)
        % Estrae i dati della connessione corrente dalla cell array.
        old_id1 = connections{i, 1}; % Nodo di partenza (ID vecchio)
        old_id2 = connections{i, 2}; % Nodo di arrivo (ID vecchio)
        type = connections{i, 3};    % Tipo geometrico ('line' o 'arc')
        extra = connections{i, 4};   % Dati aggiuntivi (es. centro dell'arco)
        L_seg = connections{i, 5};   % Lunghezza desiderata per il singolo elemento finito
        p_id = connections{i, 6};    % Proprietà fisica assegnata a questa connessione
        
        % Usa il "dizionario" lookup per trovare i nuovi ID corrispondenti.
        start_node = lookup(old_id1);
        end_node = lookup(old_id2);
        
        % Estrae le coordinate (X,Y) del nodo iniziale e finale dalla matrice riordinata.
        p1 = final_nodes(start_node, 5:6);
        p2 = final_nodes(end_node, 5:6);
        
        % Decide quale funzione di discretizzazione chiamare in base al tipo di geometria.
        if strcmp(type, 'line')
            % Se è una linea dritta, chiama la sottofunzione per interpolare lungo un segmento.
            [inter_nodes, beam_conn, next_id] = discretize_line(p1, p2, start_node, end_node, L_seg, p_id, next_id);
        else 
            % Se è un arco, estrae le coordinate del centro e il senso di percorrenza.
            center = extra(1:2);
            cw = extra(3); % Flag orario (clockwise) o antiorario
            % Chiama la sottofunzione per interpolare lungo una circonferenza.
            [inter_nodes, beam_conn, next_id] = discretize_arc(p1, p2, center, cw, start_node, end_node, L_seg, p_id, next_id);
        end
        
        % Accoda i nuovi nodi intermedi generati alla lista totale dei nodi.
        final_nodes = [final_nodes; inter_nodes];
        % Accoda le nuove travi generate alla lista totale delle travi.
        final_beams = [final_beams; beam_conn];
    end
    
    % ================================================
    % --- 3. ESTRAZIONE VETTORE ID ELEMENTI ---
    % ================================================
    % nargin = number of arguments in input. Se l'utente non ha passato 'target_prop_id'...
    if nargin < 6
        % ...restituisce un vettore colonna con tutti gli ID generati (da 1 fino al numero totale di travi).
        elements_output = (1:size(final_beams, 1))';
    else
        % ...altrimenti, cerca quali travi hanno la proprietà richiesta (terza colonna di final_beams) 
        % e restituisce gli indici di riga (che corrispondono agli ID degli elementi).
        elements_output = find(final_beams(:, 3) == target_prop_id);
    end
    
    % ====================================
    % --- 4. SCRITTURA FILE .INP ---
    % ====================================
    % Scrive l'intestazione per i Nodi. L'esclamativo (!) di solito indica un commento nel file inp.
    fprintf(fid, '*NODES\n! id c1 c2 c3 x y\n');
    for i = 1:size(final_nodes, 1)
        % %-6d = intero allineato a sinistra, %10.4f = numero con virgola (4 decimali). Stampa riga per riga.
        fprintf(fid, '%-6d %-2d %-2d %-2d  %10.4f %10.4f\n', final_nodes(i,:));
    end
    fprintf(fid, '*ENDNODES\n\n'); % Chiude il blocco nodi
    
    % Scrive l'intestazione e i dati per gli elementi Trave (Beams).
    fprintf(fid, '*BEAMS\n! id node1 node2 prop_id\n');
    for i = 1:size(final_beams, 1)
        % Qui stampa 'i' come ID della trave, seguito dai due nodi e dall'ID proprietà.
        fprintf(fid, '%-6d %-6d %-6d %-6d\n', i, final_beams(i,:));
    end
    fprintf(fid, '*ENDBEAMS\n\n');
    
    % Scrive l'intestazione e i dati per le Proprietà dei materiali/sezioni.
    fprintf(fid, '*PROPERTIES\n! id m EA EJ\n');
    for i = 1:size(properties_list, 1)
        % %12.4e stampa in notazione scientifica (es. 2.1000e+11), perfetto per EA ed EJ che sono numeri molto grandi.
        fprintf(fid, '%-4d %10.4f %12.4e %12.4e\n', properties_list(i,:));
    end
    fprintf(fid, '*ENDPROPERTIES\n');
    
    % Chiude il file salvandolo fisicamente sul disco.
    fclose(fid);
    
    % Messaggio a schermo per informare l'utente che il processo è terminato.
    fprintf('File "%s" generato: %d nodi totali, %d travi totali.\n', filename, size(final_nodes,1), size(final_beams,1));
end


% ==========================================================
% --- SOTTOFUNZIONI DI DISCRETIZZAZIONE GEOMETRICA ---
% ==========================================================

% Genera nodi ed elementi intermedi lungo un segmento rettilineo.
function [new_nodes, new_beams, next_id] = discretize_line(p1, p2, id1, id2, L_target, p_id, next_id)
    % Calcola il numero di segmenti: divide la distanza totale (norm) per la lunghezza target e arrotonda.
    % Mette un limite minimo di 1 (max(1,...)) per evitare che travi piccolissime generino 0 segmenti.
    n_seg = max(1, floor(norm(p2-p1) / L_target)+1);
    
    new_nodes = []; 
    new_beams = []; 
    curr = id1; % Imposta il nodo iniziale come nodo "corrente" per il primo elemento.
    
    for i = 1:n_seg
        if i < n_seg
            % Se non siamo all'ultimo segmento, dobbiamo creare un nuovo nodo intermedio.
            % Calcola la posizione interpolando linearmente tra p1 e p2 in base alla frazione (i/n_seg).
            pos = p1 + (i/n_seg)*(p2-p1);
            
            % Crea il nodo: imposta i vincoli a 0 (nodo libero internamente) e assegna le coordinate pos.
            new_nodes = [new_nodes; next_id, 0, 0, 0, pos];
            
            nxt = next_id;           % Il nodo di arrivo per la trave corrente è il nodo appena creato.
            next_id = next_id + 1;   % Incrementa il contatore ID globale per il prossimo nodo.
        else
            % Se siamo all'ultimo segmento, il nodo di arrivo è semplicemente il nodo finale (id2) pre-esistente.
            nxt = id2; 
        end
        % Crea l'elemento trave connettendo il nodo corrente al nodo successivo, assegnando la proprietà p_id.
        new_beams = [new_beams; curr, nxt, p_id]; 
        curr = nxt; % Il nodo di arrivo diventa il nodo di partenza per l'elemento successivo nel ciclo.
    end
end

% Genera nodi ed elementi intermedi lungo un arco di circonferenza.
function [new_nodes, new_beams, next_id] = discretize_arc(p1, p2, center, ~, id1, id2, L_target, p_id, next_id)
    % Calcola il raggio misurato rispetto al punto di inizio e al punto di fine.
    % (In un arco perfetto dovrebbero essere uguali, ma qui si calcolano entrambi per sicurezza).
    r1 = norm(p1 - center); 
    r2 = norm(p2 - center);
    
    % Calcola l'angolo (fase) in radianti del punto iniziale e finale rispetto al centro usando atan2 (copre 360°).
    t1 = atan2(p1(2)-center(2), p1(1)-center(1));
    t2 = atan2(p2(2)-center(2), p2(1)-center(1));
    
    % TRUCCO MATEMATICO: Uso della fase complessa per trovare l'angolo di apertura corretto (diff).
    % Trasformando gli angoli in numeri complessi e dividendoli (o sottraendo gli esponenti), 
    % e poi estraendo l'angolo risultante, si ottiene sempre la distanza angolare più breve tra i due vettori, 
    % aggirando il problema del "salto" matematico da +180° a -180°.
    diff = angle(exp(1i*(t2 - t1)));
    
    % Calcola quanti segmenti creare: converte l'angolo in lunghezza d'arco usando il raggio medio, divide per L_target e arrotonda.
    n_seg = max(1, floor(abs(diff) * ((r1+r2)/2) / L_target)+1);
    
    new_nodes = []; 
    new_beams = []; 
    curr = id1; % Imposta il nodo iniziale come nodo "corrente".
    
    for i = 1:n_seg
        if i < n_seg
            s = i/n_seg; % Parametro adimensionale da 0 a 1 che definisce a che punto dell'arco siamo.
            ang = t1 + s * diff; % Interpola l'angolo.
            r_interp = r1 + s * (r2 - r1); % Interpola il raggio (utile se p1 e p2 non hanno esattamente lo stesso raggio, creando una spirale).
            
            % Trasforma le coordinate polari (r_interp, ang) in coordinate cartesiane e le somma al centro dell'arco.
            pos = center + [r_interp*cos(ang), r_interp*sin(ang)];
            
            % Crea il nuovo nodo intermedio (libero, senza vincoli).
            new_nodes = [new_nodes; next_id, 0, 0, 0, pos];
            
            nxt = next_id; 
            next_id = next_id + 1;
        else
            % Per l'ultimo segmento chiude l'arco sul nodo finale noto.
            nxt = id2; 
        end
        % Definisce la topologia del piccolo elemento trave orientato sull'arco.
        new_beams = [new_beams; curr, nxt, p_id]; 
        curr = nxt;
    end
end