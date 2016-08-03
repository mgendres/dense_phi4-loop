/*
 *  mod.c
 */

int mod(int j, int k) 
{
  int ans;
  if (j<0) {
    ans=mod(j+k,k);
  } else {
    ans=j%k;
  }
  return ans;
}
