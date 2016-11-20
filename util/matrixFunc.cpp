
// Assume size of out is dealt with before this function
void matMult(mat in1, mat in2, mat out){
  if (in1.col!=in2.row) return;
  for(int k=0; k<in2.col; k++){
    for(int i=0; i<in1.row; i++){
      for(int j=0; j<in1.col; j++){
        out.matrix[i][k] += in1.matrix[i][j] * in2.matrix[j][i];
      }
    }
  }
}

// Computes transpose of matrix in
void transpose(mat in, mat out){
  float outMat[in.col][in.row];
  for(int i=0; i<in.col; i++){
    for(int j=0; j<in.row; j++){
      outMat[j][i] = in.matrix[i][j];
    }
  }
  out = mat(in.col,in.row,&outMat);
}

// Subtract
void subtractCentroid(mat in, Point center){
  int cols = in.col;
  for (int i = 0; i<cols; i++){
    for (int j=0; j<3; j++){
      in.matrix[i][j] = center[j];
    }
  }
}