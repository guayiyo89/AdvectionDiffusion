for fh in mesh.faces():
    nCell = fh.idx()
    aVt1 = []

    if mesh.is_boundary(fh):

        #guardo los vertices
        for vh in mesh.fv(fh):
            punto = mesh.point(vh)
            aVt1.append(punto)
        
        
        # gaurdo los puntos de interes
        for i in range(3):
            vecina = link_cell_to_cell[nCell][i]

            if vecina == -1:
                a = aVt1[i - 1]
                b = aVt1[i]

                #OK
                for fh in mesh.faces():
                    nCell_p = fh.idx()
                    aVt2 = []

                    if mesh.is_boundary(fh):
                    
                        for vh in mesh.fv(fh):
                            punto1 = mesh.point(vh)
                            aVt2.append(punto1)

                        for j in range(3):
                            vecina2 = link_cell_to_cell[nCell_p][j]
                            nEdge_p = link_cell_to_edge[nCell_p][j]

                            nEdge_p = int(nEdge_p)
                            vaY = aCenint_y[nEdge_p]
                            vaX = aCenint_x[nEdge_p]

                            if vecina2 == -1:
                                c = aVt2[j -1]
                                d = aVt2[j]

                                if c[0] == b[0] and d[0] == a[0]:
                                    #print('Res: ', nCell_p,'c: ', c, 'd: ',d)
                                    
                                    if vaY == nMax_Cy or vaY == nMin_Cy:                                
                                        celdaPar[nCell] = nCell_p
                                        ladoPar[nCell] = nEdge_p

                                if c[1] == b[1] and d[1] == a[1]:

                                    if vaX == nMax_Cx or vaX == nMin_Cx:                                
                                        celdaPar[nCell] = nCell_p
                                        ladoPar[nCell] = nEdge_p