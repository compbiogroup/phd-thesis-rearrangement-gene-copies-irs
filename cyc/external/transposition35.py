#!/bin/pypy3
# encoding: utf-8

import math
import sys

############################################################################
############### Do not need to sort a real permutation #####################
############################################################################

## This class was developed to sort cycle-graph configurations. We
## apply several rules in order to make the configuration closer to
## the identity. If no rule is found, so the program returns with no
## error message. 

## If the input cycle configuration is a full component, then we
## guarantee that the final permutation is the identity.

class Intergenic_Rev :
    def __init__(self, cycles, wgray, wblack) :
        self.input_cycles  = str(cycles).replace(" ", "")
        self.input_wgray   = str(wgray).replace(" ", "")[1:-1]
        self.input_wblack  = str(wblack).replace(" ", "")[1:-1]
                        
        self.graph      = cycle_configuration_graph(cycles,
                                                    wgray,
                                                    wblack)        
        self.randomized = False

    def get_num_balanced(self, vertices) :
        num_balanced = 0
        for vertice in vertices :
            if vertice[0].grays == vertice[0].blacks :
                num_balanced = num_balanced + 1
        return num_balanced
    
    def sort(self) :
        sequence = []
        graph = self.graph

        graph.calculate_cycles()
        _, vertices = graph.get_cycles(want_vertices = True)
        num_balanced = self.get_num_balanced(vertices)

        lbc = graph.n - num_balanced
        lowerb = (1.0*lbc)/2.0
        while True :
            operations  =  self.search(graph)
            if not operations :
                break

            for op in operations :
                if callable(op[0]) :
                    operations = operations + op[0](*tuple(op[1:]))
                else :
                    op = self.format_node_to_operation(op)

                    if len(op) == 3 :
                        graph.reversal(op[0],op[1],op[2])
                        sequence.append(tuple(op))
                    elif (op[0] == op[1]) or  (op[1] == op[2]):
                        graph.imove(op[0],op[1],op[2],
                                            op[3],op[4],op[5])
                        sequence.append(tuple(op))                    
                    else:
                        graph.transposition(op[0],op[1],op[2],
                                            op[3],op[4],op[5])
                        sequence.append(tuple(op))


        graph.calculate_cycles()
        _, vertices = graph.get_cycles(want_vertices = True)
        num_balanced = 0
        for vertice in vertices :
            if vertice[0].grays == vertice[0].blacks :
                num_balanced = num_balanced + 1

        if num_balanced == graph.n :
            if (((1.0*len(sequence))/lowerb) <= 4) :
                print(len(sequence))
                # print('%s %s %s %d %s %f' % (self.input_cycles,
                #                self.input_wgray,
                #                self.input_wblack,
                #                len(sequence),
                #                str(sequence).replace(" ", ""),
                #                len(sequence)/lowerb))
            else :
                sys.exit('%s %s %s %d %s %f ERROR-LOWER-BOUND-HIGHER' % (self.input_cycles,
                               self.input_wgray,
                               self.input_wblack,
                               len(sequence),
                               str(sequence).replace(" ", ""),
                               len(sequence)/lowerb))
        else :
            sys.exit('%s %s %s %d %s %f ERROR-NOT-SORTED' % (self.input_cycles,
                               self.input_wgray,
                               self.input_wblack,
                               len(sequence),
                               str(sequence).replace(" ", ""),
                               len(sequence)/lowerb))

    def format_node_to_operation(self,operation) :
        if operation :
            if len(operation) == 3 :
                i = int(round(float(operation[0].index + 1) / 2))
                j = int(round(float(operation[1].index + 1) / 2))
                return (i, j-1,operation[2])
            else :
                i = int(round(float(operation[0].index + 1) / 2))
                j = int(round(float(operation[1].index + 1) / 2))
                k = int(round(float(operation[2].index + 1) / 2))

                if i > j :
                    i,j = j,i
                if j > k :
                    j,k = k,j
                if i > j :
                    i,j = j,i
                
                return (i,j,k,int(operation[3]),int(operation[4]),int(operation[5]))



    ## Para transposições
    def search(self, graph) :
        graph.calculate_cycles()
        operations = None

        ## Uma operação de quebra.
        if not operations :
            operations = self.search_break_transposition(graph)
            if operations :
                return operations

        ## Caso simples, simplesmente procuramos um ciclo
        ## desbalanceado com mais peso nas arestas pretas e nos
        ## livramos desse peso.
        if not operations :
            operations = self.search_nontrivial_unbalanced_transposition(graph)
            if operations :         
                return [operations]

        ## Caso simples, temos um ciclo trivial com mais pesos nas
        ## arestas pretas do que nas arestas cinzas. Simplesmente
        ## jogamos esse peso extra para outro ciclo usando duas
        ## transposições.
        if not operations :
            operations = self.search_two_trivial_unbalanced_transposition(graph)
            #operations = self.search_trivial_unbalanced_imove(graph)
            if operations :          
                return operations

        ## Nesse caso, se já temos um ciclo balanceado, devemos
        ## transformar esse ciclo em orientado para poder usar o
        ## primeiro caso. O método parece complicado porque ao gerar o
        ## ciclo orientado não podemos perder ciclos balanceados.
        if not operations :
            operations = self.search_shuffling_balanced_transposition(graph)
            if operations :
                return [operations]

        ## Caso simples, temos apenas ciclos pares de tamanho 2. Achar
        ## dois cruzados e aplicar uma transposição que une os dois. O
        ## ciclo resultante será orientado e balanceado, então não
        ## preciso me preocupar em embaralhar os pesos
        if not operations :
            operations = self.search_2cycle_balanced_transposition(graph)
            if operations :
                return [operations]
            
            
    ## Para reversões
    def search_rev(self, graph) :
        graph.calculate_cycles()
        operations = None
            
        if not operations :
            ## We will merge two unbalanced cycles in order to create
            ## a big balanced cycle.
            operations = self.search_merge_balanced_reversal(graph)
            if operations :
                return self.format_node_to_operations(operations)

        if not operations :
            ## We will merge two unbalanced cycles and we guarantee
            ## that we will create on cycle that can be broken again
            ## in order to create a balanced cycle.
            operations = self.search_merge_unbalanced_reversal(graph)
            if operations :
                return self.format_node_to_operations(operations)            

        if not operations :
            ## We will merge two unbalanced cycles and we guarantee
            ## that we will create on cycle that can be broken again
            ## in order to create a balanced cycle.
            operations = self.search_non_break_reversal(graph)
            if operations :
                return self.format_node_to_operations(operations) 


##################################################
############### Several Steps ####################
##################################################


    ## Receive two nodes in the same cycle.
    def __check_convergency(self, graph, v_i, v_j) :
        if v_i.cycle == v_j.cycle :
            
            node = v_i.ap.ac 
            while node.index != v_i.index :
                node = node.ap.ac 
                if node.index == v_j.index :                    
                    return True
            return False
        sys.exit ("ERRO: Edges in different cycles")

    ## Search crossing edges. The point of this method is to find a
    ## reversal that transforms two black edges that are convergent
    ## into two black edges that are divergent. I will usually put
    ## these two black edges as parameters.
    def __search_crossing_edges(self, graph, v_i, v_j, any_cycle = False) :
        ## b is the rightmost and in the right side of the black edge
        a,b = v_i, v_j
        if a.index > b.index :
            a,b = b,a
        if a.index % 2 == 0 :
            a = a.ap
        if b.index % 2 == 0 :
            b = b.ap

        half1 = {}
        half2 = {}
        
        ## starting the search
        node = graph.end_node
        while node != b :
            half1[node.cycle] = node
            node = node.ap.ab
                    
        node = b.ap.ab
        while node != a :
            half2[node.cycle] = node
            node = node.ap.ab
        
        node = a.ap.ab
        while node :
            half1[node.cycle] = node
            node = node.ap.ab

        for key in list(half1.keys()) :
            if key in half2 :
                a,b = half1[key], half2[key]
                if ( key != v_i.cycle or
                     self.__check_convergency(graph, a, b) ) :
                    if a.index > b.index :
                        a,b = b,a
                    return a,b

    ## This is an auxiliary method called in the next method for
    ## simplicity. Remember we are talking about reversals. In order
    ## to understand this method, we advise you to draw the four cases
    ## in a piece of paper and try to code.
    def __compute_balance(self, v_i, v_j, balance, available, go_right):
        ## Check if the first cycle can be
        ## balanced. Keep in mind that balance
        ## refers to the cycle that will have v_i.
        if (0 <= balance <= available) :
            a, b = v_j, v_i
            if (v_i.index < v_j.index) :
                a,b = v_i, v_j
            exchange = balance - a.wp
            
            if go_right : ## We need balance in
                          ## the right                                    
                return (a,b, balance - b.wp) 
            else : ## We need balance in the left
                return (a, b, a.wp - balance) 
    
    ## This should be applied on black edges v_i and v_j in the same
    ## cycle. The inputs partial_gray and partial_black represents the
    ## gray and black edges that goes from v_i to v_j starting from a
    ## gray edge. In this case, it is going to be the gray edges in a
    ## cycle that would be formed if the reversal was applied
    def __get_balance(self, v_i, v_j, partial_gray, partial_black) :
                        
        ## E necessario desenhar os casos no papel para entender a
        ## linha abaixo.
        v_i_right     = (v_i.index > v_i.ap.index)
        balance      = partial_gray - partial_black

        grays  = v_i.grays  ## Sum of gray  edges
        blacks = v_i.blacks ## Sum of black edges

        available=v_i.wp + v_j.wp ## Amount of
                                  ## weight that
                                  ## can go to any
                                  ## cycle.

        go_right = v_i_right
        op = self.__compute_balance(v_i, v_j, balance, available, go_right)
        if op :
            return op
        
        go_right     = not  v_i_right                            
        remain_gray  = grays  - partial_gray 
        remain_black = blacks - partial_black - available
        balance     = remain_gray - remain_black
        op = self.__compute_balance(v_i, v_j, balance, available, go_right)
        if op :
            return op
        

    ## use this just for transpositions
    def __sum_weight_in_segment(self, init, end) :
        node = init
        blacks = 0
        grays  = 0
        while node.index != end.index :
            blacks = blacks + node.wp
            node   = node.ap
            grays  = grays  + node.wc
            node   = node.ac
        return blacks, grays



    def find_all_triples(self, cycle) :
        triples = []
        k = cycle[0]

        while True :  ## Loop over k
            i = k.ap.ac
            while True : ## Loop over i
                j = i.ap.ac                
                while True : ## Loop over j
                    if (k.index > j.index > i.index) :
                        triples.append([i, j, k] )
                    j = j.ap.ac
                    if j.index == k.index :
                        break
                i = i.ap.ac
                if i.ap.ac.index == k.index :
                    break            
            k = k.ap.ac
            if k.index == cycle[0].index :
                break
        #for triple in triples :
        #    print("Triples: ", triple[0].index, triple[1].index, triple[2].index)        
        return(triples)
    
    def find_all_triples1(self, cycle) :
        triples = []
        for k in range(0, len(cycle), 2) :
            for i in range(k+2, len(cycle), 2) :
                for j in range(i+2, len(cycle), 2) :
                    ni,nj,nk = (cycle[i], cycle[j], cycle[k])
                    if (nk.index > nj.index > ni.index) :
                        triples.append([ni, nj, nk] )
        #for triple in triples :
        #    print("Triples: ", triple[0].index, triple[1].index, triple[2].index)
        return triples


                            ## No método abaixo, vamos atribuir o peso wi
                        ## para a posição 0 da i-ésima aresta
                        ## preta. Os outros pesos são configurados em
                        ## cadeia usando gray_weights e black_weights.
    def compute_score(self, gray_ws, black_ws, ws, i, wi) :
        ws[i][0] = wi ## Início da cadeia.
        ws[i][1] = black_ws[i] - wi ## Restrição
                                    ## forte
        ws[(i+1)%3][1] = gray_ws[i] - wi ## Restrição
                                         ## fraca.                            
        ws[(i+1)%3][0] = black_ws[(i+1)%3] - ws[(i+1)%3][1]
        
        if (ws[i][0] < 0 or
            ws[i][1] < 0 or
            ws[(i+1)%3][0] < 0  or
            ws[(i+1)%3][1] < 0) :
            return -1
        
    
        ## Já colocamos todos os pesos que tinham
        ## alguma restrição. Agora temos uma
        ## aresta preta que devemos otimizar.
        best_weight, best_score = -1, 100000000000000
        for last_weight in range(0, black_ws[(i+2)%3]+1) :
            complementar = black_ws[(i+2)%3] - last_weight
            ws[(i+2)%3][0] = last_weight
            ws[(i+2)%3][1] = complementar
            score = (abs(gray_ws[(i+1)%3] - (ws[(i+1)%3][0] + ws[(i+2)%3][1])) + 
                     abs(gray_ws[(i+2)%3] - (ws[(i+2)%3][0] + ws[(i+3)%3][1])))
            if last_weight >= 0 and complementar >= 0 :
                if score < best_score :
                    best_weight = last_weight
                    best_score  = score
        ws[(i+2)%3][0] = best_weight
        ws[(i+2)%3][1] = black_ws[(i+2)%3] - ws[(i+2)%3][0]
        return best_score


    
    ## Here we deal with both balanced and unbalanced triples. We
    ## apply the first that appears.
    def search_break_transposition(self, graph) :
        _, vertices = graph.get_cycles(want_vertices = True)
        for cycle in vertices :
            if (cycle[0].direction == 1 and ## Check if it is convergent
                cycle[0].size      >= 3 and ## Check if it is long
                cycle[0].blacks    >= cycle[0].grays ## Check if it is
                                                     ## breakable,
                                                     ## which means we
                                                     ## need more
                                                     ## weight in
                                                     ## black edges
                                                     ## than gray
                                                     ## edges
            ) :
                triples = self.find_all_triples(cycle)
                #for tr in triples:
                #    print tr[0].index, tr[1].index, tr[2].index
              
                for transp in triples :
                    ## Let us compute the segment weights between
                    ## the edges that will be broken. We will call
                    ## it gray_weights, even though it may not
                    ## represent weights
                    gray_weights = []
                    ## Pesos das arestas pretas
                    black_weights = []
                    ## Here I will place the transp weights. A
                    ## soma de cada bloco tem que dar o peso da
                    ## aresta preta correspondente.
                    transp_weights = [[0,0],[0,0],[0,0]] 

                    ## Computing the gray black weights.
                    for i in range(len(transp)) :
                        init, end = transp[i], transp[(i+1)%3]
                        blacks, grays = self.__sum_weight_in_segment(init, end)
                        blacks = blacks - init.wp
                        # transp_weights[i][0] = init.wp
                        black_weights.append(init.wp)
                        gray_weights.append(grays - blacks)
 
                    ## Não entendi o código anterior, vou criar um
                    ## que faz mais sentido.
                    
                    ## Passo 1: Devo escolher uma aresta preta com
                    ## peso maior que ou igual à aresta cinza
                    ## (função g(.) na verdade) incidente a ela
                    ## pela esquerda. Podemos provar que existe
                    ## pelo menos uma satisfazendo esse critério
                    ## por contradição.

                    # Trocando
                    # heaviest_black = None
                    # for i in range(3) :
                    #     if black_weights[i] >= gray_weights[i] :
                    #         heaviest_black = i
                    #         break

                    best_score  = 100000000000000
                    best_weight = 100000000000000
                    best_black  = -1
                    
                    for heaviest_black in (0,1,2) :                            
                        ## Passo 2: A aresta cinza à esquerda de
                        ## heaviest_black deve garantir um ciclo
                        ## balanceado porque a aresta preta tem peso
                        ## para dar a ela. O que devemos fazer agora,
                        ## é garantir que as duas outras arestas
                        ## pretas vão se aproximar do peso das suas
                        ## aresta cinzas. Para isso, vou criar um
                        ## score para uma atribuição de pesos. Esse
                        ## score é simplesmente a soma absoluta da
                        ## diferença de pesos entre as arestas pretas
                        ## e cinzas na configuração resultante.
    
    
                        for weight in range(0, black_weights[heaviest_black]+1) :
                            score = self.compute_score(gray_weights, black_weights,
                                                       transp_weights, heaviest_black,
                                                       weight)
    
                            if score >= 0 and score < best_score :
                                best_score  = score
                                best_weight = weight
                                best_black  = heaviest_black
                                                                            
                    if best_black != -1 and best_score == 0 :
                        self.compute_score(gray_weights, black_weights, transp_weights,
                                           best_black, best_weight)
                            
                                    
                        outcome = [
                            gray_weights[0] - (transp_weights[0][0] + transp_weights[1][1]),
                            gray_weights[1] - (transp_weights[1][0] + transp_weights[2][1]),
                            gray_weights[2] - (transp_weights[2][0] + transp_weights[0][1])
                        ]
                        
                        end_point = 0
                        for el in outcome :
                            if el == 0:
                                end_point = end_point + 1
                        
                        start_point = int(cycle[0].blacks  == cycle[0].grays)
                                                
                        #if True :
                        if end_point - start_point > 0 :
                            return [(transp[0], transp[1],
                                     transp[2],
                                    transp_weights[0][0], transp_weights[1][0],
                                     transp_weights[2][0])]
                        else :
                            return [
                                (transp[0], transp[1], transp[2],
                                 transp_weights[0][0], transp_weights[1][0], transp_weights[2][0]),
                                (self.__search_two_trivial_unbalanced_transposition,
                                #(self.__search_trivial_unbalanced_imove,
                                 transp[0], transp[1], transp[2])
                                ]

                for transp in triples :
                    ## Let us compute the segment weights between
                    ## the edges that will be broken. We will call
                    ## it gray_weights, even though it may not
                    ## represent weights
                    gray_weights = []
                    ## Pesos das arestas pretas
                    black_weights = []
                    ## Here I will place the transp weights. A
                    ## soma de cada bloco tem que dar o peso da
                    ## aresta preta correspondente.
                    transp_weights = [[0,0],[0,0],[0,0]] 

                    ## Computing the gray black weights.
                    for i in range(len(transp)) :
                        init, end = transp[i], transp[(i+1)%3]
                        blacks, grays = self.__sum_weight_in_segment(init, end)
                        blacks = blacks - init.wp
                        # transp_weights[i][0] = init.wp
                        black_weights.append(init.wp)
                        gray_weights.append(grays - blacks)

                    ## Não entendi o código anterior, vou criar um
                    ## que faz mais sentido.
                    
                    ## Passo 1: Devo escolher uma aresta preta com
                    ## peso maior que ou igual à aresta cinza
                    ## (função g(.) na verdade) incidente a ela
                    ## pela esquerda. Podemos provar que existe
                    ## pelo menos uma satisfazendo esse critério
                    ## por contradição.

                    # Trocando
                    # heaviest_black = None
                    # for i in range(3) :
                    #     if black_weights[i] >= gray_weights[i] :
                    #         heaviest_black = i
                    #         break

                    best_score  = 100000000000000
                    best_weight = 100000000000000
                    best_black  = -1
                    
                    for heaviest_black in (0,1,2) :                            
                        ## Passo 2: A aresta cinza à esquerda de
                        ## heaviest_black deve garantir um ciclo
                        ## balanceado porque a aresta preta tem peso
                        ## para dar a ela. O que devemos fazer agora,
                        ## é garantir que as duas outras arestas
                        ## pretas vão se aproximar do peso das suas
                        ## aresta cinzas. Para isso, vou criar um
                        ## score para uma atribuição de pesos. Esse
                        ## score é simplesmente a soma absoluta da
                        ## diferença de pesos entre as arestas pretas
                        ## e cinzas na configuração resultante.
    
    
                        for weight in range(0, black_weights[heaviest_black]+1) :
                            score = self.compute_score(gray_weights, black_weights,
                                                       transp_weights, heaviest_black,
                                                       weight)
    
                            if score >= 0 and score < best_score :
                                best_score  = score
                                best_weight = weight
                                best_black  = heaviest_black
                                                                            
                    if best_black != -1 :
                        self.compute_score(gray_weights, black_weights, transp_weights,
                                           best_black, best_weight)                        
                                    
                        outcome = [
                            gray_weights[0] - (transp_weights[0][0] + transp_weights[1][1]),
                            gray_weights[1] - (transp_weights[1][0] + transp_weights[2][1]),
                            gray_weights[2] - (transp_weights[2][0] + transp_weights[0][1])
                        ]
                        
                        end_point = 0
                        for el in outcome :
                            if el == 0:
                                end_point = end_point + 1
                        
                        start_point = int(cycle[0].blacks  == cycle[0].grays)
                        
                        #if True :
                        if end_point - start_point > 0 :
                            return [(transp[0], transp[1],
                                     transp[2],
                                    transp_weights[0][0], transp_weights[1][0],
                                     transp_weights[2][0])]
                        else :
                            return [
                                (transp[0], transp[1], transp[2],
                                 transp_weights[0][0], transp_weights[1][0], transp_weights[2][0]),
                                (self.__search_two_trivial_unbalanced_transposition,
                                #(self.__search_trivial_unbalanced_imove,
                                 transp[0], transp[1], transp[2])
                                ]
                    #print("ERROR: Avisar o Ulisses Imediatamente " +  graph.to_string())

    ## Here I get a nonoriented cycle and I get edges in the open
    ## gates. As a second criterium, I select only the heaviest edges.
    def __search_interleaving_edges(self, graph, v_i, v_j, v_k) :
        ## c is the rightmost and in the right side of the black
        ## edge. I first put everybody in its place.
        a,b,c = v_i, v_j, v_k
        if a.index > b.index :
            a,b = b,a
        if b.index > c.index :
            b,c = c,b
        if a.index > b.index :
            a,b = b,a
            
        if a.index % 2 == 0 :
            a = a.ap
        if b.index % 2 == 0 :
            b = b.ap
        if c.index % 2 == 0 :
            c = b.ap
            
        ## Now I create the halves
            
        half0 = {} ## outside
        half1 = {}
        half2 = {}

        ## starting the search
        graph.clean_visit()
        node = graph.end_node
        while node != c :
            node.visit = 0
            if node.cycle in half0 and node.wp < half0[node.cycle].wp :
                pass
            else :
                half0[node.cycle] = node
            node = node.ap.ab
                
        ## starting the search
        node = c.ap.ab
        while node != b :
            node.visit = 1
            if node.cycle in half1 and node.wp < half1[node.cycle].wp :
                pass
            else :
                half1[node.cycle] = node
            node = node.ap.ab
                    
        node = b.ap.ab
        while node != a :
            node.visit = 2            
            if node.cycle in half2 and node.wp < half2[node.cycle].wp :
                pass
            else :
                half2[node.cycle] = node
            node = node.ap.ab
        
        node = a.ap.ab
        while node :
            node.visit = 0            
            if node.cycle in half0 and node.wp < half0[node.cycle].wp :
                pass
            else :                
                half0[node.cycle] = node
            node = node.ap.ab

        ## Let us first try interleave
        for key in list(half0.keys()) :
            if key in half1 and key in half2 :
                a,b,c = half0[key], half2[key], half1[key]
                if ( key != v_i.cycle ) :
                    if a.index > c.index :
                        return b, c, a
                    else :
                        return a, b, c

        ## Let us clean the dictionary, there might be many unitary
        ## cycles or cycles in only one half. We do not need them.
        for half in (half0, half1, half2) :
            for key in list(half.keys()) :
                num = 0
                if key in half0 :
                    num = num + 1
                if key in half1 :
                    num = num + 1
                if key in half2 :
                    num = num + 1
                ## It is in only one half. So we remove it
                if num == 1 :
                    half.pop(key, None)
                    
        ## Let us try shattered and balanced. A balanced is in only
        ## two blocks, so this loop could be improved, however I will
        ## let it as it is for now.
        for half in (half0, half1, half2) :
            for key in list(half.keys()) :
                if (key != v_i.cycle and half[key].grays == half[key].blacks) :
                    ## Now I need to find two edges in different regions
                    ## such that sum of black edges is bigger than the
                    ## balance in both sides.
                    keep_searching = True
                    found          = False
                    init = half[key]
                    end  = init.ap.ac
                    while keep_searching :
                        while end != init and keep_searching :
                            if end.visit != init.visit :
                                blacks, grays = self.__sum_weight_in_segment(init, end)
                                blacks = blacks - init.wp
                                if 0 <= grays - blacks <= init.wp + end.wp :
                                    keep_searching = False
                                    found          = True
                                    break
                            end = end.ap.ac
                        if keep_searching :
                            init = init.ap.ac
                            if init == half[key] :
                                keep_searching = False
                                found          = False
        
                    if found :
                        a,b,c = init, end, None
                        visited = (init.visit, end.visit)
                        if 0 not in visited :
                            c = half0.popitem()[1]
                        elif 1 not in visited :
                            c = half1.popitem()[1]
                        else :
                            c = half2.popitem()[1]
                            
                        # In this case we return out of order
                        # if a.index > b.index :
                        #     a,b = b,a
                        # if b.index > c.index :
                        #     b,c = c,b
                        # if a.index > b.index :
                        #     a,b = b,a                        
                        
                        return a, b, c
                            
        ## Let us try shattered. There is no balanced here.
        for key in list(half0.keys()) :
            if (key != v_i.cycle) :
                a,b,c = None, None, None
                if key in half1 :
                    a,c = half0[key], half1[key]
                    b   = half2.popitem()[1] ## Primeiro node do dicionário
                if key in half2 :
                    a,b = half0[key], half2[key]
                    c   = half1.popitem()[1] ## Primeiro node do dicionário
                if ( a and b and c ) :
                    if a.index > c.index :
                        return b, c, a
                    else :
                        return a, b, c
        sys.exit("""ERROR: this line should never be reached, are you sure you removed the oriented cycles? """)

    ## Here we have a transposition where a and b are in the same
    ## cycle and need to get rid of weight to send to c. We need to
    ## guaratee that the cycle that will be generated by C=(...a...b)
    ## and will not contain c is balanced.
    def __compute_weight_to_send_to_other_cycle(self, a, b, c) :
        if a.index > b.index :
            a,b = b,a
            
        ## when c is between a and b, the inner
        ## edge of the blackmost cycle will be
        ## alone in a cycle. So, I will replace a
        ## and b again.
        if a.index < c.index < b.index :
            blacks, grays = self.__sum_weight_in_segment(b, a)            
            blacks = blacks - b.wp
            balance = grays - blacks
            if 0 <= balance <= b.wp + a.wp :                            
                wa = min(balance, a.wp)
                wb = balance - wa
                wa = a.wp - wa
                return a, c, b, wa, c.wp, wb
            else :
                return None
        else :
            blacks, grays = self.__sum_weight_in_segment(a, b)
            blacks = blacks - a.wp
            balance = grays - blacks
            if 0 <= balance <= a.wp + b.wp :                            
                wa = min(balance, a.wp)
                wb = b.wp - (balance - wa)
                if c.index < a.index :
                    return (c, a, b, c.wp, wa, wb)
                else :
                    return (a, b, c, wa, wb, c.wp)                
            else :
                return None
        
    ## At this point, we know that no oriented balanced cycle exists.
    def search_shuffling_balanced_transposition(self, graph) :
        _, vertices = graph.get_cycles(want_vertices = True)
        for cycle in vertices :
            if (cycle[0].direction == 1 and ## Check if it is convergent
                cycle[0].size      >=  3 and ## Check if it is long
                cycle[0].blacks    == cycle[0].grays 
            ) :
                ## The following link works because no oriented
                ## balanced exists.
                transp = self.__search_interleaving_edges(
                    graph,
                    cycle[4],
                    cycle[2],
                    cycle[0])

                ## If it is a shuffling transposition, we can do it.
                if transp[0].cycle == transp[1].cycle == transp[2].cycle :
                    return (transp[0], transp[1], transp[2],
                            transp[0].wp,
                            transp[1].wp,
                            transp[2].wp)
                ## If the cycles are unbalanced, we can do it.
                if ( (transp[0].blacks != transp[0].grays) and
                     (transp[1].blacks != transp[1].grays) and
                     (transp[2].blacks != transp[2].grays) ) :
                    return (transp[0], transp[1], transp[2],
                            transp[0].wp,
                            transp[1].wp,
                            transp[2].wp)
                
                ## At least one of them is balanced. So we need to be
                ## careful at this point. We need to move a balanced
                ## amount of weight.  We know that
                ## self.__search_interleaving_edges places the
                ## balanced cycle first, which is great. Also, the
                ## returned edges are already the ones we can break,
                ## we just need to compute weights.
                a,   b,  c = transp[0], transp[1], transp[2]
                return self.__compute_weight_to_send_to_other_cycle(a, b, c)

    
            
    ## At this point, we know that no oriented balanced cycle exists.
    def search_nontrivial_unbalanced_transposition(self, graph) :
        _, vertices = graph.get_cycles(want_vertices = True)
        for cycle in vertices :
            if (cycle[0].direction == 1 and ## Check if it is convergent
                cycle[0].size      >= 2 and ## Check if it is nontrivial
                cycle[0].blacks    > cycle[0].grays 
            ) :
                ## Let us first find another cycle where we will put
                ## the extra weight.
                c = None
                for cycle1 in vertices :
                    if (cycle1[0].blacks  < cycle1[0].grays) :
                        c = cycle1[0]
                        break
                    
                ## Let us find a long cycle
                init = cycle[0]
                end  = init.ap.ac
                while True : ## I will be locked in this infinity loop
                             ## to depure if something fails
                    while end != init :
                        op = self.__compute_weight_to_send_to_other_cycle(init, end, c) 
                        if op :
                            return op
                        end = end.ap.ac
                    init = init.ap.ac
                    end  = init.ap.ac
                    if init == cycle[0] :
                        break


    def __search_trivial_unbalanced_transposition(self, node0, node1, node2) :
        a,b,c = node0, node1, node2
        
        #ANDRE changed the following two lines
        wa, wb, wc = 0, 0, c.wp
        #wa, wb, wc = a.wc, b.wp, c.wp

        if a.index > b.index :
            a,b = b,a
            wa, wb = wb, wa
        if b.index > c.index :
            b,c = c,b
            wb, wc = wc, wb
        if a.index > b.index :
            a,b = b,a
            wa, wb = wb, wa

        #BEGIN ANDRE ADDED THIS HELPS TO GET BETTER RESULTS BUT IT IS NOT IN THE PROOF
        ## If there are three cycles such that the sum of grays and blacks are equal,
        ## we are going to change A, B and C to balance the three with 2 operations
        if (a.wp + b.wp + c.wp) == (a.wc + b.wc + c.wc) and a.blacks == a.wp and b.blacks == b.wp and c.blacks == c.wp:
            if node0 == a:
                wa, wb, wc = min(a.wp,max(a.wc,a.wc+b.wc-b.wp)), b.wp, 0
                awpnew = wa      + 0
                bwpnew = a.wp-wa + 0
                cwpnew = b.wp    + c.wp
                wa2, wb2, wc2 = a.wc, bwpnew, cwpnew+bwpnew-c.wc
            elif node0 == b:
                wa, wb, wc = min(a.wp,a.wc), b.wp-min(b.wp,max(0,a.wc-a.wp)), c.wp-min(c.wp,c.wc)
                awpnew = wa      + b.wp-wb
                cwpnew = c.wp-wc + wb
                bwpnew = a.wp-wa + wc
                wa2, wb2, wc2 = min(a.wc,awpnew), bwpnew-(a.wc-awpnew), b.wc-(awpnew-min(a.wc,awpnew))
            else:
                wa, wb, wc = a.wp, 0, c.wp-min(c.wp,c.wc+b.wc)
                awpnew = a.wp    + b.wp
                cwpnew = c.wp-wc  + 0
                bwpnew = wc + 0
                wa2, wb2, wc2 = min(awpnew,a.wc), 0, cwpnew-c.wc
        ## Otherwise, just move positive balance from cycle0 (a)
        ## to node1 (b) and leave node2 (c) unchanged
        else:
        #if True:
            if node2 == a:
                wa, wb, wc = a.wp, b.wp, 0
                if node0 == b:
                    wa2, wb2, wc2 = a.wp, 0, b.wc
                else:
                    wa2, wb2, wc2 = a.wp, 0, b.wp+c.wp-c.wc
            elif node2 == b:
                wa, wb, wc = 0, b.wp, c.wp
                if node0 == a:
                    wa2, wb2, wc2 = 0, a.wp+c.wp-a.wc, b.wp
                else:
                    wa2, wb2, wc2 = 0, c.wc, b.wp
            else:
                wa, wb, wc = a.wp, 0, 0
                if node0 == a:
                    wa2, wb2, wc2 = a.wc, 0, 0 #fix
                else:
                    wa2, wb2, wc2 = a.wp+b.wp-b.wc, 0, 0 #fix
                
                #wa, wb, wc = a.wp, b.wp, 0
        #END ANDRE ADDED THIS

        return [[a, b, c, wa, wb, wc],[b, a, c, wa2, wb2, wc2]]
        

    ## At this point, we know that no blackheaviest cycle is bigger
    ## than 1. Therefore we have little ammunition, our last resource
    ## is to send the extra weight somewhere.
    def search_trivial_unbalanced_transposition(self, graph) :
        _, vertices = graph.get_cycles(want_vertices = True)
        for cycle in vertices :
            if (cycle[0].blacks    > cycle[0].grays) :

                ## Search for other gray edges where we will send the
                ## extra weight.
                node1, node2 = None, None
                node = graph.end_node                
                while node :
                    #####ANDRE changed != to <
                    if (node.blacks < node.grays and
                        node.cycle  != cycle[0].cycle) :
                        if   not node1 :
                            node1 = node
                        else :
                            node2 = node
                            break
                    node = node.ap.ab

                ## There is only two unbalanced cycles and they are
                ## trivial, so I pick any cycle as node2.
                node = graph.end_node
                while not node2 :
                    if (node.index != node1.index and
                        node.index != cycle[0].index) :
                        node2 = node
                    node = node.ap.ab

                return self.__search_trivial_unbalanced_transposition(cycle[0], node1, node2)


    def __search_two_trivial_unbalanced_transposition(self, node0, node1, node2) :  
        a,b,c = node0, node1, node2

        #ANDRE changed the following two lines
        wa, wb, wc = 0, 0, c.wp
        #wa, wb, wc = a.wc, b.wp, c.wp

        if a.index > b.index :
            a,b = b,a
            wa, wb = wb, wa
        if b.index > c.index :
            b,c = c,b
            wb, wc = wc, wb
        if a.index > b.index :
            a,b = b,a
            wa, wb = wb, wa

        ### Verifica se node1 tambem eh negativo, se sim apenas node2 eh positivo
        ### e manda todo o excesso pra ele
        if node1.blacks > node1.grays:
            if node2 == a:
                wa = a.wp
                wb = b.wc
                wc = c.wp - c.wc
                wa2 = a.wp
                wb2 = 0
                wc2 = b.wc
                return [[a, b, c, wa, wb, wc],[b, a, c, wa2, wb2, wc2]]

            elif node2 == b:
                wa = a.wp
                wb = b.wp
                wc = 0
                wa2 = a.wc
                wb2 = 0
                wc2 = b.wp+c.wp-c.wc
                return [[a, b, c, wa, wb, wc],[b, a, c, wa2, wb2, wc2]]

            else:
                wa = a.wc
                wb = b.wp
                wc = 0
                wa2 = a.wc
                wb2 = a.wp-a.wc
                wc2 = b.wc
                return [[a, b, c, wa, wb, wc],[b, a, c, wa2, wb2, wc2]]


        ### Se nao, verifca se node2 eh do mesmo ciclo que node1 ou eh balanceado
        ### entao apenas joga de node0 para node1 deixando node2 intacto
        elif (node2.cycle == node1.cycle) or (node2.grays == node2.blacks):
            if node2 == a:
                if node0 == b: #node1 == c
                    wa = a.wp
                    wb = b.wp
                    wc = 0
                    wa2 = a.wp
                    wb2 = 0
                    wc2 = b.wc
                    return [[a, b, c, wa, wb, wc],[b, a, c, wa2, wb2, wc2]]

                else: #node0 == c node1 == b
                    wa = a.wp
                    wb = b.wp
                    wc = 0
                    wa2 = a.wp
                    wb2 = 0
                    wc2 = wb+b.grays-b.blacks
                    return [[a, b, c, wa, wb, wc],[b, a, c, wa2, wb2, wc2]]

            elif node2 == b:
                if node0 == a: #node1 == c
                    wa = a.wc
                    wb = b.wp
                    wc = 0
                    wa2 = a.wc
                    wb2 = a.wp-wa
                    wc2 = b.wp
                    return [[a, b, c, wa, wb, wc],[b, a, c, wa2, wb2, wc2]]

                else: #node0 == c node1 == a
                    wa = a.wp
                    wb = 0
                    wc = c.wp - c.wc
                    wa2 = a.wp
                    wb2 = 0
                    wc2 = 0
                    return [[a, b, c, wa, wb, wc],[b, a, c, wa2, wb2, wc2]]

            else: #node2 == c
                if node0 == a: #node1 == b
                    wa = a.wp
                    wb = 0
                    wc = 0
                    wa2 = a.wc
                    wb2 = 0
                    wc2 = 0
                    return [[a, b, c, wa, wb, wc],[b, a, c, wa2, wb2, wc2]]

                else: #node0 == b node1 == a
                    wa = a.wp
                    wb = 0
                    wc = 0
                    wa2 = a.wp+b.wp-b.wc
                    wb2 = 0
                    wc2 = 0
                    return [[a, b, c, wa, wb, wc],[b, a, c, wa2, wb2, wc2]]

        ### Se node2 nao eh negativo nem balanceado, entao ele eh positivo.
        ### Como node2 nao esta no mesmo ciclo que node1, precisamos balancear
        ### node0 e node1, jogar o resto em node2
        else:
            if node2.cycle != node1.cycle:
                if node0 == a:
                    wa = a.wc + (b.grays - b.blacks)
                    wb = 0
                    wc = 0
                    wa2 = a.wc
                    wb2 = a.wp - wa
                    wc2 = 0
                    return [[a, b, c, wa, wb, wc],[b, a, c, wa2, wb2, wc2]]

                elif node0 == b:
                    wa = a.wp
                    wb = c.grays-c.blacks
                    wc = 0
                    wa2 = b.wc
                    wb2 = 0
                    wc2 = 0
                    return [[a, b, c, wa, wb, wc],[b, a, c, wa2, wb2, wc2]]

                else:
                    wa = a.wp
                    wb = b.wp
                    wc = c.wp - c.wc - (b.grays-b.blacks)
                    wa2 = a.wp
                    wb2 = 0
                    wc2 = b.grays-b.blacks
                    return [[a, b, c, wa, wb, wc],[b, a, c, wa2, wb2, wc2]]

        sys.exit("ERRO Ops, nenhuma operacao encontrada...")
        sys.exit()

    ## At this point, we know that no blackheaviest cycle is bigger
    ## than 1. Therefore we have little ammunition, our last resource
    ## is to send the extra weight somewhere. The trick here is that
    ## we can find either two negative, or two positives or one of each
    ## and create two new balanced cycles
    def search_two_trivial_unbalanced_transposition(self, graph) :
        _, vertices = graph.get_cycles(want_vertices = True)
        for cycle in vertices :
            if (cycle[0].blacks    > cycle[0].grays) :

                ## Search another cycle that is also negative to be node1
                node1, node2 = None, None
                node = graph.end_node                
                while node :
                    if (node.blacks > node.grays and node.cycle != cycle[0].cycle) :
                        if   not node1 :
                            node1 = node
                    node = node.ap.ab

                ## If node1 is None, then let us choose a positive cycle to be node1
                if not node1:
                    node = graph.end_node                
                    while node :
                        if (node.blacks < node.grays and node.cycle != cycle[0].cycle) :
                            if not node1 :
                                node1 = node
                        node = node.ap.ab
                
                ## Let us find now a positive cycle (that is different from node1) to be node2
                node = graph.end_node                
                while node :
                    if (node.blacks < node.grays and node.cycle != node1.cycle) :
                        if   not node2 :
                            node2 = node
                    node = node.ap.ab

                ## There is only two unbalanced cycles, so I pick any balanced cycle as node2.
                node = graph.end_node
                while not node2 :
                    if (node.index != node1.index and
                        node.index != cycle[0].index) :
                        node2 = node
                    node = node.ap.ab

                return self.__search_two_trivial_unbalanced_transposition(cycle[0], node1, node2)

    def __search_trivial_unbalanced_imove(self, node0, node1, node2) :
        a,b,c = node0, node1, node2

        #imove na direita
        if a.index > b.index :
            return [[b, a, a, b.wp, a.wc, a.wp],]
        #imove na esquerda
        else :
            return [[a, a, b, a.wc, a.wp, b.wp],]

    ## At this point, we know that no blackheaviest cycle is bigger
    ## than 1. Therefore we have little ammunition, our last resource
    ## is to send the extra weight somewhere.
    def search_trivial_unbalanced_imove(self, graph) :
        _, vertices = graph.get_cycles(want_vertices = True)
        for cycle in vertices :
            if (cycle[0].blacks    > cycle[0].grays) :

                ## Search for other gray edges where we will send the
                ## extra weight.
                node1, node2 = None, None
                node = graph.end_node                
                while node :
                    #####ANDRE changed != to <
                    if (node.blacks < node.grays and
                        node.cycle  != cycle[0].cycle) :
                        if   not node1 :
                            node1 = node
                        else :
                            node2 = node
                            break
                    node = node.ap.ab

                ## There is only two unbalanced cycles and they are
                ## trivial, so I pick any cycle as node2.
                node = graph.end_node
                while not node2 :
                    if (node.index != node1.index and
                        node.index != cycle[0].index) :
                        node2 = node
                    node = node.ap.ab

                return self.__search_trivial_unbalanced_imove(cycle[0], node1, node2)
                
    ## Search crossing edges. The point of this method is to find a
    ## reversal that transforms two black edges that are convergent
    ## into two black edges that are divergent. I will usually put
    ## these two black edges as parameters.
    def __search_crossing_edges_for_transposition(self, graph, v_i, v_j) :
        ## b is the rightmost and in the right side of the black edge
        a,b = v_i, v_j
        if a.index > b.index :
            a,b = b,a
        if a.index % 2 == 0 :
            a = a.ap
        if b.index % 2 == 0 :
            b = b.ap

        half1 = {}
        half2 = {}
        
        ## starting the search
        node = graph.end_node
        while node != b :
            half1[node.cycle] = node
            node = node.ap.ab
                    
        node = b.ap.ab
        while node != a :
            half2[node.cycle] = node
            node = node.ap.ab        
        node = a.ap.ab
        while node :
            half1[node.cycle] = node
            node = node.ap.ab

        for key in list(half1.keys()) :
            if key in half2 :
                a,b = half1[key], half2[key]
                if ( key != v_i.cycle  ) :
                    if a.index > b.index :
                        a,b = b,a
                    return a,b

            
    def search_2cycle_balanced_transposition(self, graph) :
        _, vertices = graph.get_cycles(want_vertices = True)
        for cycle in vertices :
            if (cycle[0].size   ==  2) :
                a,b = cycle[0].ap.ac, cycle[0]
                c,d = self.__search_crossing_edges_for_transposition(graph, a, b)
                op = self.__compute_weight_to_send_to_other_cycle(a, b, c) 
                if op :
                    return op
                
    
    ## Procuro um ciclo divergente e aplico imediatamente uma reversao
    ## em duas de suas arestas pretas (divergentes) que balanceiam
    ## pelo menos um dos ciclos.
    def search_break_reversal(self, graph) :
        _, vertices = graph.get_cycles(want_vertices = True)
        for cycle in vertices :
            ## Check if it is divergent
            if cycle[0].direction == 2 :
                ## Look for two divergent black edges that can be
                ## broken. This can be probably replaced by the
                ## biggest black edge.
                for i in range(1, len(cycle) - 2, 2) :
                    v_i           = cycle[i]
                    
                    partial_gray  =  v_i.wc
                    partial_black =  0
                    for j in range(i+2, len(cycle), 2)  :
                        v_j           = cycle[j]

                        if (v_i.index%2!=v_j.index%2): ## Check divergence
                            op = self.__get_balance(v_i, v_j, partial_gray, partial_black)
                            if op :
                                return op
                            
                        ## Preparing next iteration    
                        partial_gray  = partial_gray  + v_j.wc
                        partial_black = partial_black + v_j.wp


    def search_non_break_reversal(self, graph) :
        _, vertices = graph.get_cycles(want_vertices = True)
        for cycle in vertices :
                ## Look for two divergent black edges that can be
                ## broken. This can be probably replaced by the
                ## biggest black edge.
                for i in range(1, len(cycle) - 2, 2) :
                    v_i           = cycle[i]
                    
                    partial_gray  =  v_i.wc
                    partial_black =  0
                    for j in range(i+2, len(cycle), 2)  :
                        v_j           = cycle[j]

                        op = self.__get_balance(v_i, v_j, partial_gray, partial_black)
                        if op :
                            return op
                            
                        ## Preparing next iteration    
                        partial_gray  = partial_gray  + v_j.wc
                        partial_black = partial_black + v_j.wp


                        
    def search_merge_balanced_reversal(self, graph) :
        _, vertices = graph.get_cycles(want_vertices = True)
        unbalanced   = []
        for cycle in vertices :
            if cycle[0].grays != cycle[0].blacks :
                unbalanced.append(cycle[0])
        for i in range(len(unbalanced)) :
            cycle_i = unbalanced[i]
            for j in range(i+1, len(unbalanced)) :
                cycle_j = unbalanced[j]
                if (cycle_i.grays  + cycle_j.grays ==
                    cycle_i.blacks + cycle_j.blacks) :
                    a,b = cycle_i, cycle_j
                    if (a.index > b.index) :
                        a,b = b,a
                    return a, b, 0

    def search_crossing_reversal(self, graph) :
        _, vertices = graph.get_cycles(want_vertices = True)
        for cycle in vertices :
            for i in range(1, len(cycle) - 2, 2) :
                v_i           = cycle[i]
                partial_gray  =  v_i.wc
                partial_black =  0
                for j in range(i+2, len(cycle), 2)  :
                    v_j           = cycle[j]
                    op = self.__get_balance(v_i, v_j, partial_gray, partial_black)

                    if op :
                        op = self.__search_crossing_edges(graph, v_i, v_j)
                        if op :
                            return op[0], op[1], 0
                            
                    ## Preparing next iteration    
                    partial_gray  = partial_gray  + v_j.wc
                    partial_black = partial_black + v_j.wp

        

    def search_merge_unbalanced_reversal(self, graph) :
        _, vertices = graph.get_cycles(want_vertices = True)
        unbalanced   = []
        for cycle in vertices :
            if cycle[0].grays != cycle[0].blacks :

                biggest = cycle[0]
                for i in range(0, len(cycle), 2)  :
                    if (cycle[i].wp > biggest.wp) :
                        biggest = cycle[i]
                unbalanced.append( cycle[i]  )
            
                
        for i in range(len(unbalanced)) :
            cycle_i = unbalanced[i]
            for j in range(i+1, len(unbalanced)) :
                cycle_j = unbalanced[j]
                
                if (
                        ## remover pretas
                        cycle_i.blacks - cycle_i.wp <= cycle_i.grays <= cycle_i.blacks or
                        cycle_j.blacks - cycle_j.wp <= cycle_j.grays <= cycle_j.blacks or
                        ## inserir blacks do outro
                        cycle_i.blacks <= cycle_i.grays <= cycle_i.blacks + cycle_j.wp or
                        cycle_j.blacks <= cycle_j.grays <= cycle_j.blacks + cycle_i.wp) :
                    a,b = cycle_i, cycle_j
                    if (a.index > b.index) :
                        a,b = b,a
                    return a, b, 0
        
    
############################################################################
############ Do not need to represent a real permutation ###################
############################################################################

## This class represents an unoriented cycle-graph. Keep in mind that
## this class is not the same as used for the sorting by transposition
## problem. Here, the white edges are not the reverse of the black
## edges as before. That is because the black edges can go from righ
## to left as well as to left to right.

## For instance, assume a permutation (5, 2, ...). The configuration
## in the graph would be.

## node(0).ap  = -5      ## node(0).ab  = Null
## node(-5).ab = 5       ## node(-5).ap = 0  
## node(5).ap  = -2      ## node(5).ab  = -5 
## node(-2).ab = 2       ## node(-2).ap = 5

class cycle_configuration_graph() :
    def __init__(self, cycles,
                 weight_gray ,
                 weight_black) :

        if (sum(weight_gray) != sum(weight_black)) :
            sys.exit("ERRO Different intergenic regions")
            sys.exit()


        ## The number of cycles was never calculated
        self.__num_cycles           = False
        self.__num_odd_cycles       = False
        self.__first_indice_shift   = 0 # Tell us the index of the
                                        # edge that is pointed by
                                        # begin_node. If it is
                                        # negative, than we know that
                                        # a mirror have occurred.

        ## self.n is the number of black edges. Remember this graph
        ## might not be a permutation
        self.n = 0
        for cycle in cycles :
            self.n = self.n + len(cycle)
        n = self.n
          
        # Creating nodes
        node_list = []
        node_list = [cycle_graph_node(i, False) for i in range(2*n)]
        self.begin_node = node_list[0 ]
        self.end_node   = node_list[-1]

        # Creating ap
        for i in range(0,2*n,2) :
            node_list[i  ].ap  = node_list[i+1]
            node_list[i+1].ap  = node_list[i  ]
            node_list[i  ].wp  = weight_black[int(i/2)]
            node_list[i+1].wp  = weight_black[int(i/2)]

        # Creating ab
        for i in range(1,2*n-1,2) :
            node_list[i  ].ab  = node_list[i+1]
            node_list[i+1].ab  = node_list[i  ]

        # Creating ac 
        for cycle in cycles :
            for i in range(len(cycle)) :               
                front, back  = -1,-1 # from left of the black edge.
                j = (i+1)%len(cycle)

                if cycle[i] > 0 :
                    front = 2*( cycle[i]  -1 )
                else :
                    front = 2*(-cycle[i]) -1

                if cycle[j] > 0 :
                    back = 2*( cycle[j]) -1
                else :
                    back = 2*(-cycle[j]  -1 )

                node_list[front].ac = node_list[back]
                node_list[back].ac  = node_list[front]

        self.set_values()
        for i in range(0,2*n,1) :
            positive =  min(abs(node_list[i].value),
                            abs(node_list[i].ac.value))
            
            weight = weight_gray[positive]
            node_list[i].wc     = weight
            node_list[i].ac.wc  = weight

############################################################                
################ Rearrangement Operations  #################
############################################################                
    ## 
    def transposition(self, i, j, k, weight_i, weight_j, weight_k) :
        node_i = None
        node_j = None
        node_k = None

        unp_i, unp_j, unp_k = 0,0,0

        count     = 0
        unp_count = 0

        node = self.begin_node
        # Find the black edges
        while node :
            count = count + 1
            if not node.padded :
                unp_count = unp_count + 1

            if count == i :
                node_i = node
                unp_i  = unp_count
            if count == j :
                node_j = node
                unp_j  = unp_count
            if count == k :
                node_k = node
                unp_k  = unp_count
            node = node.ap.ab
        
        if (weight_i <  0 or weight_i > node_i.wp) :
            sys.exit("ERRO: Peso inconsistente no lado esquerdo: i (max = %s: w = %s)" % (
                node_i.wp,
                weight_i
            ) )
            sys.exit()

        if (weight_j < 0 or weight_j > node_j.wp) :
            sys.exit("ERRO: Peso inconsistente no meio:  j (max = %s: w = %s)" % (
                node_j.wp,
                weight_j
            ) )
            sys.exit()

        if (weight_k < 0 or weight_k > node_k.wp) :
            sys.exit("ERRO: Peso inconsistente no lado direito:  k (max = %s: w = %s)" % (
                node_k.wp,
                weight_k
            ) )
            sys.exit()


        ## Weights that will be exchanged.
        lefti, righti = weight_i, node_i.wp - weight_i
        leftj, rightj = weight_j, node_j.wp - weight_j
        leftk, rightk = weight_k, node_k.wp - weight_k
            
        # Change the edges
        node_i.ap.ap = node_k
        node_j.ap.ap = node_i
        node_k.ap.ap = node_j

        
        aux       = node_i.ap
        node_i.ap = node_j.ap
        node_j.ap = node_k.ap
        node_k.ap = aux

        # Updating weights: Node_i keeps its place, Node_j moved to
        # the third black edge and Node_k is now in the middle.
        node_i.wp = node_i.ap.wp = lefti  + rightj
        node_k.wp = node_k.ap.wp = righti + leftk
        node_j.wp = node_j.ap.wp = rightk + leftj
        
        self.__num_cycles           = False
        self.__num_odd_cycles       = False
        self.reset_indices()
        return unp_i, unp_j, unp_k
        
    def reversal(self, i, j, weight) :
        node_i = None
        node_j = None

        unp_i, unp_j = 0,0

        count     = 0
        unp_count = 0
        ## These variables are used to return the transposition
        ## applied as integer. ret_i, and ret_j will be the same
        ## as i,j if we do not have ignored black edges
        node = self.begin_node
        # Find the black edges

        while node :
            count = count + 1
            if not node.padded :
                unp_count = unp_count + 1

            if count == i :
                node_i = node
                unp_i  = unp_count
            if count == j+1 :
                node_j = node
                unp_j  = unp_count - 1
            node = node.ap.ab


        if (weight >  0 and weight > node_i.wp) :
            sys.exit("ERRO: Peso inconsistente no lado esquerdo: i (max = %s: w = %s)" % (
                node_i.wp,
                weight
            ) )
            sys.exit()

        if (weight < 0 and abs(weight) > node_j.wp) :
            sys.exit("ERRO: Peso inconsistente no lado direito:  j (max = %s: w = %s)" % (
                node_j.wp,
                weight
            ) )
            sys.exit()

        nweight_i = node_i.wp - weight
        nweight_j = node_j.wp + weight            
            
        # Change the edges
        node_i.ap.ap = node_j.ap
        node_i.ap.wp = nweight_j        
        node_j.ap.ap = node_i.ap
        node_j.ap.wp = nweight_j        

        node_i.ap  = node_j
        node_i.wp  = nweight_i        
        node_j.ap  = node_i
        node_j.wp  = nweight_i

        self.__num_cycles           = False
        self.__num_odd_cycles       = False
        self.__num_balanced_cycles  = False                
        self.reset_indices()
        return unp_i, unp_j

    def imove(self, i, j, k, weight_i, weight_j, weight_k) :
        node_i = None
        node_j = None
        node_k = None

        unp_i, unp_j, unp_k = 0,0,0

        count     = 0
        unp_count = 0

        node = self.begin_node
        # Find the black edges
        while node :
            count = count + 1
            if not node.padded :
                unp_count = unp_count + 1

            if count == i :
                node_i = node
                unp_i  = unp_count
            if count == j :
                node_j = node
                unp_j  = unp_count
            if count == k :
                node_k = node
                unp_k  = unp_count
            node = node.ap.ab

        # Check if arguments are consistent
        if (i == j and i != k) :
            if (weight_i <  0 or weight_i > node_i.wp) :
                sys.exit("ERRO: Peso inconsistente no start point do lado esquerdo do intergenic move: i (max = %s: w = %s)" % (
                    node_i.wp,
                    weight_i
                ) )
                sys.exit()
            elif (weight_j <=  weight_i or weight_j > node_i.wp) :
                sys.exit("ERRO: Peso inconsistente no end point do lado esquerdo do intergenic move: j (max = %s: w = %s)" % (
                    node_i.wp,
                    weight_j
                ) )
                sys.exit()

        elif (i != j and j == k) :
            if (weight_j <  0 or weight_j > node_k.wp) :
                sys.exit("ERRO: Peso inconsistente no start point do lado direito do intergenic move: i (max = %s: w = %s)" % (
                    node_k.wp,
                    weight_j
                ) )
                sys.exit()
        elif (weight_k <=  weight_j or weight_k > node_k.wp) :
                sys.exit("ERRO: Peso inconsistente no end point do lado direito do intergenic move: j (max = %s: w = %s)" % (
                    node_k.wp,
                    weight_k
                ) )
                sys.exit()

        else :
            sys.exit("ERRO: Intergenic move com indices i = j = k")
            sys.exit()

        # Updating weights:
        if (i == j) :
            wimove = weight_j - weight_i
            node_i.wp = node_i.wp - wimove
            node_k.wp = node_k.wp + wimove
        else :
            wimove = weight_k - weight_j
            node_k.wp = node_k.wp - wimove
            node_i.wp = node_i.wp + wimove
        
        self.__num_cycles           = False
        self.__num_odd_cycles       = False
        self.reset_indices()
        return unp_i, unp_j, unp_k
        

############################################################                
######################## Bounds  ###########################
############################################################                
    def upper_bound(self) :
        return self.n - self.num_cycles()

    def lower_bound(self) :
        return int(math.ceil( (self.n - self.num_odd_cycles())/2.0 ))

############################################################                
################### Auxiliary Methods  #####################
############################################################                
    def print_indices(self) :
        node = self.begin_node

        while node :
            print("%i, padded = %s" % (node.index, node.padded))
            node = node.ap
            print("%i, padded = %s" % (node.index, node.padded))
            node = node.ab

    def to_string(self) :
        self.clean_visit()        

        node = self.end_node        

        str_cycles = "("

        while node :
            if not node.visit :
                cycle_node = node
                cycle      = "("

                while not cycle_node.visit :
                    wc = cycle_node.wc
                    wp = cycle_node.wp
                    if cycle_node.index % 2 == 0 :
                        cycle = cycle + "%i<%i,%i>," % (-(cycle_node.index+2)/2,wp,wc)
                    else :
                        cycle = cycle + "%i<%i,%i>," % (+(cycle_node.index+1)/2,wp,wc )


                    cycle_node.visit = True
                    cycle_node = cycle_node.ap
                    cycle_node.visit  = True
                    cycle_node = cycle_node.ac
                    
                cycle = cycle[:-1] + "),"
                str_cycles = str_cycles + cycle
            node = node.ap.ab
        str_cycles = str_cycles[:-1] + ")"
        return str_cycles

            

    def get_cycles(self, want_vertices = False) :
        self.clean_visit()        

        node = self.end_node        
        cycles    = []
        vertices  = []

        while node :
            if not node.visit :
                cycle_node  = node
                cycle       = []
                cycle_nodes = []

                while not cycle_node.visit :
                    if cycle_node.index % 2 == 0 :
                        cycle.append( -(cycle_node.index+2)/2 )
                    else :
                        cycle.append( +(cycle_node.index+1)/2 )


                    cycle_node.visit = True
                    cycle_nodes.append(cycle_node)
                    cycle_node = cycle_node.ap
                    cycle_node.visit  = True
                    cycle_nodes.append(cycle_node)
                    cycle_node = cycle_node.ac
                    
                cycles.append(tuple(cycle))
                vertices.append(cycle_nodes)

            node = node.ap.ab
        if want_vertices :
            return tuple(cycles), vertices
        return tuple(cycles)

    ## Here we try to give values to the nodes to check if they are a
    ## permutation.
    def set_values(self) :
        node = self.begin_node        

        node.value = 0
        node = node.ac
        node.value = -1

        while node.ab :
            node = node.ab
            node.value = -(node.ab.value)
            node = node.ac

            if node.ac.value > 0 :
                node.value = -(node.ac.value + 1)
            else :
                node.value = -(node.ac.value - 1)
    
    def reset_indices(self) :
        node  = self.begin_node 
        count = 0

        while node :
            node.index = count
            node       = node.ap
            
            count      = count + 1
            
            node.index = count
            node       = node.ab
            
            count      = count + 1

    def num_cycles(self) :
        if type(self.__num_cycles) == bool :
            self.calculate_cycles()
        return self.__num_cycles

    def num_odd_cycles(self) :
        if type(self.__num_odd_cycles) == bool :
            self.calculate_cycles()
        return self.__num_odd_cycles


    def calculate_cycles(self) :
        cycles, vertices = self.get_cycles(want_vertices = True)
        num_cycles = len(cycles)
        num_odd = 0

        for cycle in cycles :
            if len(cycle) % 2 == 1 :
                num_odd = num_odd + 1

        self.__num_cycles           = num_cycles
        self.__num_odd_cycles       = num_odd

        for i in range(len(cycles)) :            
            cycle = cycles[i]
            size = len(cycle)

            direction = 1
            for el in cycle :
                if (el < 0) :
                    direction = 2
                    break

            blacks   = 0
            grays    = 0
            vertice_set = vertices[i]
            for vertex in vertice_set :
                grays   = grays   + vertex.wc
                blacks  = blacks  + vertex.wp
            grays  = grays  / 2
            blacks = blacks / 2

            vertice_set = vertices[i]
            for vertex in vertice_set :
                vertex.size      = size
                vertex.grays     = grays
                vertex.blacks    = blacks
                vertex.direction = direction                
                #vertex.kind  = kind
                vertex.cycle = i

    def clean_visit(self) :
        node = self.begin_node
        
        while node :
            node.visit = False
            node       = node.ap
            node.visit = False
            node       = node.ab

        
#####################################################################
################## REPRESENTS A NODE OF A GRAPH #####################
#####################################################################

class cycle_graph_node :
    #value  : stores pi_i    
    #index  : stores the black edge i, 0 <= i <= n+1
    #ap     : points to the record that stores pi_{i-1}, 1 <= i <= n+1
    #ab     : points to the record that stores pi_{i+1}, 0 <= i <= n
    #ac     : points to the record that stores i + 1,    0 <= i <= n
    #rac    : reverse of the ac
    #visit  : indicates if the record has been visited
    #cend   : stores the black edge i_k in the canonical representation
    def __init__(self, index, padded) :
        self.index, self.value    = index, 0        
        self.padded   = padded
        self.cycle    = 0
        ## Kind:
        ## 0 = unset
        ## 1 = small cycle             (1)
        ## 2 = oriented    convergent  (3, 1, 2)
        ## 3 = even        convergent  (2,  1)
        ## 4 = even        divergent   (2, -1)
        ## 5 = nonoriented convergent  (3, 2, 1)
        ## 6 = oriented    divergent   (3,-2,-1), (3,1,-2), (3,-1, 2)
        ## 7 = nonoriened  divergent   (3,2,-1),  (3,-2,1), (3,-1,-2)
        self.kind     = 0

        ## 0 = unset 
        ## i = num of black edges
        self.size        = 0 

        ## 0   = unset or balanced
        ## +-i = (gray - black) weights
        self.blacks      = 0
        self.grays       = 0
        
        ## 0 = unset
        ## 1 = convergent
        ## 2 = divergent
        self.direction = 0

        
        self.ap, self.ab, self.ac = 0,0,0
        self.visit  = False
        self.wc                   = 0
        self.wp                   = 0
        
#        self.cend = 0

def get_position(permutation) :
    n = len(permutation)-2
    position    = [-1 for i in range(0, n+2)]
    for i in range(0, n+2) :
        position[abs(permutation[i])] = i
    return position

def get_rightmost_element(cycle, position) :
    max_position = 0
    for i in range(len(cycle)) :
        if position[cycle[i]] > position[cycle[max_position]] :
            max_position = i
    return max_position

## The unordered cycle starts with a gray edge, we order them by
## making it start with the rightmost black edge.
def order_cycle(cycle, position) :
    index = get_rightmost_element(cycle, position)
    new   = []
    new.append(cycle[index])

    if index % 2 == 0 :
        iter_el  = (index-1) % len(cycle)
        while iter_el != index :
            new.append(cycle[iter_el])
            iter_el = (iter_el-1) % len(cycle)
    else :
        iter_el  = (index+1) % len(cycle)
        while iter_el != index :
            new.append(cycle[iter_el])
            iter_el = (iter_el+1) % len(cycle)
    return new

def canonical_representation(cycle, position) :
    cycle     = order_cycle(cycle, position)
    canonical = []

    for i in range(0,len(cycle),2) :
        if position[cycle[i]] < position[cycle[i+1]] :
            black = -position[cycle[i+1]]
            canonical.append(black )
        else :
            black = position[cycle[i]]
            canonical.append(black)
    return canonical


def construct_str_cycle(permutation) :
    n = len(permutation)-2

    # permutation = [0] + permutation + [n+1]
    position    = [-1 for i in range(0, n+2)]
    #sign        = [-1 for i in range(0, n+2)]


    for i in range(1, n+2) :
        position[abs(permutation[i])] = i
    #    sign    [abs(permutation[i])] = permutation[i] / abs(permutation[i])

    ## 1 if the gray edge i,-(i+1) was never used.
    gray_available     = [1 for i in range(0, n+1)]
    #black_available    = [1 for i in range(0, n+1)]

    cycles = []

    for i in range(0, n+1) :

        ## 
      if gray_available[i] :

        start     = i
        cycle = [start]

        end   = start
        positive  = True
        
        while True :
            
            ## Will be used later, it says if after walking through
            ## the black edge we are in the right or in the left
            is_vertex_left = None

            if positive :
                ## Gray edge: we are looking for the edge ( end,-(end+1) )
                gray_available[end] = gray_available[end] - 1
                end = end + 1
                cycle.append(end)

                ## Black edge: we are at the vertex -end.
                if permutation[position[end]] > 0 :
                    # If the sign in that position is positive, than
                    # -end is in the left (cycle graph)                    
                    end = abs(permutation[position[end]-1])
                    is_vertex_left = False
                
                else :
                    # If the sign in that position is negative, than
                    # -end is in the right (cycle graph)
                    end = abs(permutation[position[end]+1])
                    is_vertex_left = True
            else :
                ## Gray edge: we are looking for the edge ( -end, end-1  )
                end = end - 1                                 ##  Note we swapped
                gray_available[end] = gray_available[end] - 1 ##  these lines
                cycle.append(end)


                ## Black edge: we are at the vertex +end.
                if permutation[position[end]] > 0 :
                    # If the sign in that position is positive, than
                    # +end is in the right (cycle graph)
                    end = abs(permutation[position[end]+1])
                    is_vertex_left = True
                else : 
                    # If the sign in that position is negative, than
                    # +end is in the left (cycle graph)
                    end = abs(permutation[position[end]-1])
                    is_vertex_left = False
                    
            if end == start :
                break
            else :
                cycle.append(end)
                
                if is_vertex_left :
                    if permutation[position[end]] < 0 :
                        positive = True
                    else :
                        positive = False
                else :                    
                    if permutation[position[end]] < 0 :
                        positive = False
                    else :
                        positive = True
        cycles.append(cycle)

    int_position = get_position(permutation)
    canonicals = []

    for cycle in cycles :
        canonicals.append(canonical_representation(cycle, int_position))

    str_canonicals = str(canonicals)
    str_canonicals = str_canonicals.replace(" ", "")
    str_canonicals = str_canonicals.replace("[", "(")
    str_canonicals = str_canonicals.replace("]", ",)")
    return str_canonicals

        
if __name__ == '__main__':
    permutation = eval("[%s]" % sys.argv[1])
    str_cycles = construct_str_cycle(permutation)

    config = eval(str_cycles)
    wgray  = eval("[%s]" % sys.argv[2])
    wblack = eval("[%s]" % sys.argv[3])
    sort = Intergenic_Rev(config, wgray, wblack)
    sort.sort()

