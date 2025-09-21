#include "eispack.h"

double r8_abs(double x) {
  return x >= 0.0 ? x : -x;
}

double r8_max(double x, double y) {
  return y < x ? x : y;
}

double r8_min(double x, double y) {
  return y < x ? y : x;
}

double r8_sign(double x) {
  return x < 0.0 ? -1.0 : 1.0;
}

double r8_epsilon() {
  static double value = 2.220446049250313E-016;
  return value;
}

double pythag(double a, double b) {
  double p = r8_max(r8_abs(a), r8_abs(b));
  if (p != 0.0) {
    double r = r8_min(r8_abs(a), r8_abs(b)) / p;
    r = r * r;
    for (;;) {
      double t = 4.0 + r;
      if (t == 4.0) break;
      double s = r / t;
      double u = 1.0 + 2.0 * s;
      p = u * p;
      r = (s / u) * (s / u) * r;
    }
  }
  return p;
}

void tred1(int n, double a[], double d[], double e[], double e2[]) {
  double f, g, h, scale;
  int i, j, k, l;

  for (j = 0; j < n; j++) d[j] = a[n - 1 + j * n];
  for (i = 0; i < n; i++) a[n - 1 + i * n] = a[i + i * n];

  for (i = n - 1; 0 <= i; i--) {
    l = i - 1;
    h = 0.0;
    scale = 0.0;
    for (k = 0; k <= l; k++) scale += r8_abs(d[k]);

    if (scale == 0.0) {
      for (j = 0; j <= l; j++) {
        d[j] = a[l + j * n];
        a[l + j * n] = a[i + j * n];
        a[i + j * n] = 0.0;
      }
      e[i] = 0.0;
      e2[i] = 0.0;
      continue;
    }

    for (k = 0; k <= l; k++) d[k] /= scale;
    for (k = 0; k <= l; k++) h += d[k] * d[k];

    e2[i] = h * scale * scale;
    f = d[l];
    g = -sqrt(h) * r8_sign(f);
    e[i] = scale * g;
    h -= f * g;
    d[l] = f - g;

    if (0 <= l) {
      for (k = 0; k <= l; k++) e[k] = 0.0;

      for (j = 0; j <= l; j++) {
        f = d[j];
        g = e[j] + a[j + j * n] * f;
        for (k = j + 1; k <= l; k++) {
          g += a[k + j * n] * d[k];
          e[k] += a[k + j * n] * f;
        }
        e[j] = g;
      }

      f = 0.0;
      for (j = 0; j <= l; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      h = f / (h + h);
      for (j = 0; j <= l; j++) e[j] -= h * d[j];

      for (j = 0; j <= l; j++) {
        f = d[j];
        g = e[j];
        for (k = j; k <= l; k++) a[k + j * n] -= f * e[k] - (-g) * d[k];
      }
    }

    for (j = 0; j <= l; j++) {
      f = d[j];
      d[j] = a[l + j * n];
      a[l + j * n] = a[i + j * n];
      a[i + j * n] = f * scale;
    }
  }
}

void tred2(int n, double a[], double d[], double e[], double z[]) {
  double f, g, h, hh, scale;
  int i, j, k, l;

  for (j = 0; j < n; j++)
    for (i = j; i < n; i++) z[i + j * n] = a[i + j * n];

  for (j = 0; j < n; j++) d[j] = a[n - 1 + j * n];

  for (i = n - 1; 1 <= i; i--) {
    l = i - 1;
    h = 0.0;
    scale = 0.0;
    for (k = 0; k <= l; k++) scale += r8_abs(d[k]);

    if (scale == 0.0) {
      e[i] = d[l];
      for (j = 0; j <= l; j++) {
        d[j] = z[l + j * n];
        z[i + j * n] = 0.0;
        z[j + i * n] = 0.0;
      }
      d[i] = 0.0;
      continue;
    }

    for (k = 0; k <= l; k++) d[k] /= scale;

    h = 0.0;
    for (k = 0; k <= l; k++) h += d[k] * d[k];

    f = d[l];
    g = -sqrt(h) * r8_sign(f);
    e[i] = scale * g;
    h -= f * g;
    d[l] = f - g;

    for (k = 0; k <= l; k++) e[k] = 0.0;

    for (j = 0; j <= l; j++) {
      f = d[j];
      z[j + i * n] = f;
      g = e[j] + z[j + j * n] * f;
      for (k = j + 1; k <= l; k++) {
        g += z[k + j * n] * d[k];
        e[k] += z[k + j * n] * f;
      }
      e[j] = g;
    }

    for (k = 0; k <= l; k++) e[k] /= h;
    f = 0.0;
    for (k = 0; k <= l; k++) f += e[k] * d[k];
    hh = 0.5 * f / h;
    for (k = 0; k <= l; k++) e[k] -= hh * d[k];

    for (j = 0; j <= l; j++) {
      f = d[j];
      g = e[j];
      for (k = j; k <= l; k++) z[k + j * n] -= f * e[k] + g * d[k];
      d[j] = z[l + j * n];
      z[i + j * n] = 0.0;
    }
    d[i] = h;
  }

  for (i = 1; i < n; i++) {
    l = i - 1;
    z[n - 1 + l * n] = z[l + l * n];
    z[l + l * n] = 1.0;
    h = d[i];
    if (h != 0.0) {
      for (k = 0; k <= l; k++) d[k] = z[k + i * n] / h;
      for (j = 0; j <= l; j++) {
        g = 0.0;
        for (k = 0; k <= l; k++) g += z[k + i * n] * z[k + j * n];
        for (k = 0; k <= l; k++) z[k + j * n] -= g * d[k];
      }
    }
    for (k = 0; k <= l; k++) z[k + i * n] = 0.0;
  }

  for (j = 0; j < n; j++) d[j] = z[n - 1 + j * n];
  for (j = 0; j < n - 1; j++) z[n - 1 + j * n] = 0.0;
  z[n - 1 + (n - 1) * n] = 1.0;
  e[0] = 0.0;
}

int tql2(int n, double d[], double e[], double z[]) {
  double c, c2, c3 = 0, dl1, el1, f, g, h, p, r, s, s2 = 0, t, tst1, tst2;
  int i, ierr = 0, ii, j, k, l, l1, l2, m, mml;

  if (n == 1) return ierr;

  for (i = 1; i < n; i++) e[i - 1] = e[i];
  f = 0.0;
  tst1 = 0.0;
  e[n - 1] = 0.0;

  for (l = 0; l < n; l++) {
    j = 0;
    h = r8_abs(d[l]) + r8_abs(e[l]);
    tst1 = r8_max(tst1, h);

    for (m = l; m < n; m++) {
      tst2 = tst1 + r8_abs(e[m]);
      if (tst2 == tst1) break;
    }

    if (m != l) {
      for (;;) {
        if (30 <= j) { ierr = l + 1; return ierr; }
        j++;
        l1 = l + 1;
        l2 = l1 + 1;
        g = d[l];
        p = (d[l1] - g) / (2.0 * e[l]);
        r = pythag(p, 1.0);
        d[l]  = e[l] / (p + r8_sign(p) * r8_abs(r));
        d[l1] = e[l] * (p + r8_sign(p) * r8_abs(r));
        dl1 = d[l1];
        h = g - d[l];
        for (i = l2; i < n; i++) d[i] -= h;
        f += h;

        p = d[m];
        c = 1.0;
        c2 = c;
        el1 = e[l1];
        s = 0.0;
        mml = m - l;

        for (ii = 1; ii <= mml; ii++) {
          c3 = c2; c2 = c; s2 = s;
          i = m - ii;
          g = c * e[i];
          h = c * p;
          r = pythag(p, e[i]);
          e[i + 1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i + 1] = h + s * (c * g + s * d[i]);

          for (k = 0; k < n; k++) {
            h = z[k + (i + 1) * n];
            z[k + (i + 1) * n] = s * z[k + i * n] + c * h;
            z[k + i * n] = c * z[k + i * n] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
        tst2 = tst1 + r8_abs(e[l]);
        if (tst2 <= tst1) break;
      }
    }
    d[l] += f;
  }

  for (ii = 1; ii < n; ii++) {
    i = ii - 1;
    k = i;
    p = d[i];
    for (j = ii; j < n; j++) if (d[j] < p) { k = j; p = d[j]; }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j < n; j++) {
        t = z[j + i * n];
        z[j + i * n] = z[j + k * n];
        z[j + k * n] = t;
      }
    }
  }
  return ierr;
}

int tqlrat(int n, double d[], double e2[]) {
  double b = 0, c = 0, f, g, h, p, r, s, t;
  int i, ierr = 0, ii, j, l, l1, m, mml;

  if (n == 1) return ierr;

  for (i = 1; i < n; i++) e2[i - 1] = e2[i];
  f = 0.0;
  t = 0.0;
  e2[n - 1] = 0.0;

  for (l = 0; l < n; l++) {
    j = 0;
    h = r8_abs(d[l]) + sqrt(e2[l]);
    if (t <= h) { t = h; b = r8_abs(t) * r8_epsilon(); c = b * b; }

    for (m = l; m < n; m++) if (e2[m] <= c) break;

    if (m != l) {
      for (;;) {
        if (30 <= j) { ierr = l + 1; return ierr; }
        j++;
        l1 = l + 1;
        s = sqrt(e2[l]);
        g = d[l];
        p = (d[l1] - g) / (2.0 * s);
        r = pythag(p, 1.0);
        d[l] = s / (p + r8_abs(r) * r8_sign(p));
        h = g - d[l];
        for (i = l1; i < n; i++) d[i] -= h;
        f += h;

        g = d[m]; if (g == 0.0) g = b;
        h = g; s = 0.0; mml = m - l;

        for (ii = 1; ii <= mml; ii++) {
          i = m - ii;
          p = g * h;
          r = p + e2[i];
          e2[i + 1] = s * r;
          s = e2[i] / r;
          d[i + 1] = h + s * (h + d[i]);
          g = d[i] - e2[i] / g;
          if (g == 0.0) g = b;
          h = g * p / r;
        }
        e2[l] = s * g;
        d[l] = h;

        if (h == 0.0) break;
        if (r8_abs(e2[l]) <= r8_abs(c / h)) break;
        e2[l] = h * e2[l];
        if (e2[l] == 0.0) break;
      }
    }

    p = d[l] + f;
    for (i = l; 0 <= i; i--) {
      if (i == 0) { d[i] = p; break; }
      else if (d[i - 1] <= p) { d[i] = p; break; }
      d[i] = d[i - 1];
    }
  }
  return ierr;
}

int rs(int n, double a[], double w[], int matz, double z[]) {
  double *fv1, *fv2;
  int ierr;

  if (matz == 0) {
    fv1 = new double[n];
    fv2 = new double[n];
    tred1(n, a, w, fv1, fv2);
    ierr = tqlrat(n, w, fv2);
    delete[] fv1;
    delete[] fv2;
  } else {
    fv1 = new double[n];
    tred2(n, a, w, fv1, z);
    ierr = tql2(n, w, fv1, z);
    delete[] fv1;
  }
  return ierr;
}