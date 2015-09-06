// DC-GL-Matrix

glMatrixArrayType = typeof Float32Array != "undefined" ? Float32Array : typeof WebGLFloatArray != "undefined" ? WebGLFloatArray:Array;

var vec2 = {};

vec2.create = function(a) {
	var b = new glMatrixArrayType(2);
	if (a) {
		b[0] = a[0];
		b[1] = a[1];
	}
	return b;
}

vec2.set = function(a, b) {
	b[0] = a[0];
	b[1] = a[1];
	return b;
}

vec2.add = function(a, b, c) {
	if (!c || a == c) {
		a[0] += b[0];
		a[1] += b[1];
		return a;
	}
	c[0] = a[0] + b[0];
	c[1] = a[1] + b[1];
	return c;
}

vec2.subtract = function(a, b, c) {
	if (!c || a == c) {
		a[0] -= b[0];
		a[1] -= b[1];
		return a;
	}
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	return c;
}

vec2.negate = function(a, b) {
	b || (b = a);
	b[0] = -a[0];
	b[1] = -a[1];
	return b;
};

vec2.scale = function(a, b, c) {
	if(!c || a == c) {
		a[0] *= b;
		a[1] *= b;
		return a;
	}
	c[0] = a[0] * b;
	c[1] = a[1] * b;
	return c;
}

vec2.normalize = function(a, b) {
	b || (b = a);
	var c = a[0], d = a[1], g = Math.sqrt(c * c + d * d);
	if(g) {
		if(g == 1) {
			b[0] = c;
			b[1] = d;
			return b;
		}
	}
	else {
		b[0] = 0;
		b[1] = 0;
		return b;
	}
	g = 1 / g;
	b[0] = c * g;
	b[1] = d * g;
	return b;
}

/*
vec2.cross = function(a, b, c) {
	c || (c = a);
	var d = a[0], e = a[1];
	a = a[2];
	var g = b[0], f = b[1];
	b = b[2];
	c[0] = e * b - a * f;
	c[1] = a * g - d * b;
	c[2] = d * f - e * g;
	return c;
}
*/

vec2.length = function(a) {
	var b = a[0];
	a = a[1];
	return Math.sqrt(b * b + a * a);
}

vec2.dot = function(a, b) {
	return a[0] * b[0] + a[1] * b[1];
}

vec2.direction = function(a, b, c) {
	c || (c = a);
	var d = a[0] - b[0];
	a = a[1] - b[1];
	b = Math.sqrt(d * d + a * a);
	if (!b) {
		c[0] = 0;
		c[1] = 0;
		return c;
	}
	b = 1 / b;
	c[0] = d * b;
	c[1] = a * b;
	return c;
}

vec2.lerp = function(a, b, c, d) {
	d || (d = a);
	d[0] = a[0] + c * (b[0] - a[0]);
	d[1] = a[1] + c * (b[1] - a[1]);
	return d;
}

vec2.str = function(a) {
	return "[" + a[0] + ", " + a[1] + "]";
}

var mat3 = {};

// Original
mat3.create = function(a) {
	var b = new glMatrixArrayType(9);
	if (a) {
		b[0] = a[0];
		b[1] = a[1];
		b[2] = a[2];
		b[3] = a[3];
		b[4] = a[4];
		b[5] = a[5];
		b[6] = a[6];
		b[7] = a[7];
		b[8] = a[8];
		b[9] = a[9];
	}
	return b;
}

// Original
mat3.set = function(a, b) {
	b[0] = a[0];
	b[1] = a[1];
	b[2] = a[2];
	b[3] = a[3];
	b[4] = a[4];
	b[5] = a[5];
	b[6] = a[6];
	b[7] = a[7];
	b[8] = a[8];
	return b;
}

// Original
mat3.identity = function(a) {
	a[0] = 1;
	a[1] = 0;
	a[2] = 0;
	a[3] = 0;
	a[4] = 1;
	a[5] = 0;
	a[6] = 0;
	a[7] = 0;
	a[8] = 1;
	return a;
}

mat3.transpose = function(a,b) {
	if (!b || a == b) {
		var c = a[1], d = a[2], e = a[5];
		a[1]=a[3];
		a[2]=a[6];
		a[3]=c;
		a[5]=a[7];
		a[6]=d;
		a[7]=e;
		return a;
	}
	b[0]=a[0];
	b[1]=a[3];
	b[2]=a[6];
	b[3]=a[1];
	b[4]=a[4];
	b[5]=a[7];
	b[6]=a[2];
	b[7]=a[5];
	b[8]=a[8];
	return b;
}

// New
mat2.toMat3 = function(a, b) {
	b || (b = mat3.create());
	b[0] = a[0];
	b[1] = a[1];
	b[2] = 0;
	b[3] = a[2];
	b[4] = a[3];
	b[5] = 0;
	b[6] = 0;
	b[7] = 0;
	b[8] = 1;
	return b;
}

// Original
mat3.str = function(a) {
	return "[" + a[0] + ", " + a[1] + ", " + a[2] + ", " + a[3] + ", " + a[4] + ", " + a[5] + ", " + a[6] + ", " + a[7] + ", " + a[8] + "]";
}

// New









var mat2 = {};

mat2.create = function(a) {
var b = new glMatrixArrayType(4);
	if(a){
		b[0]=a[0];
		b[1]=a[1];
		b[2]=a[2];
		b[3]=a[3];
	}
	return b;
}

mat2.set = function(a, b) {
	b[0] = a[0];
	b[1] = a[1];
	b[2] = a[2];
	b[3] = a[3];
	return b;
}

mat2.identity = function(a) {
	a[0] = 1;
	a[1] = 0;
	a[2] = 0;
	a[3] = 1;
	return a;
}

mat2.transpose = function(a, b) {
	if(!b || a == b) {
		var c = a[1];
		a[1] = a[2];
		a[2] = c;
		return a;
	}
	b[0] = a[0];
	b[1] = a[2];
	b[2] = a[1];
	b[3] = a[3];
	return b;
}

mat2.toMat3 = function(a, b) {
	b || (b = mat3.create());
	b[0] = a[0];
	b[1] = a[1];
	b[2] = 0;
	b[3] = a[2];
	b[4] = a[3];
	b[5] = 0;
	b[6] = 0;
	b[7] = 0;
	b[8] = 1;
	return b;
}

mat2.toMat4 = function(a, b) {
	b || (b = mat4.create());
	b[0] = a[0];
	b[1] = a[1];
	b[2] = 0;
	b[3] = 0;
	b[4] = a[2];
	b[5] = a[3];
	b[6] = 0;
	b[7] = 0;
	b[8] = 0;
	b[9] = 0;
	b[10] = 1;
	b[11] = 0;
	b[12] = 0;
	b[13] = 0;
	b[14] = 0;
	b[15] = 1;
	return b;
}

mat2.str = function(a) {
	return "[" + a[0] + ", " + a[1] + ", " + a[2] + ", " + a[3] + "]";
}

mat2.determinant = function(a) {
var b=a[0],c=a[1],d=a[2],e=a[3],g=a[4],f=a[5],h=a[6],i=a[7],j=a[8],k=a[9],l=a[10],o=a[11],m=a[12],n=a[13],p=a[14];a=a[15];return m*k*h*e-j*n*h*e-m*f*l*e+g*n*l*e+j*f*p*e-g*k*p*e-m*k*d*i+j*n*d*i+m*c*l*i-b*n*l*i-j*c*p*i+b*k*p*i+m*f*d*o-g*n*d*o-m*c*h*o+b*n*h*o+g*c*p*o-b*f*p*o-j*f*d*a+g*k*d*a+j*c*h*a-b*k*h*a-g*c*l*a+b*f*l*a;
}

mat2.inverse=function(a,b){b||(b=a);var c=a[0],d=a[1],e=a[2],g=a[3],f=a[4],h=a[5],i=a[6],j=a[7],k=a[8],l=a[9],o=a[10],m=a[11],n=a[12],p=a[13],r=a[14],s=a[15],A=c*h-d*f,B=c*i-e*f,t=c*j-g*f,u=d*i-e*h,v=d*j-g*h,w=e*j-g*i,x=k*p-l*n,y=k*r-o*n,z=k*s-m*n,C=l*r-o*p,D=l*s-m*p,E=o*s-m*r,q=1/(A*E-B*D+t*C+u*z-v*y+w*x);b[0]=(h*E-i*D+j*C)*q;b[1]=(-d*E+e*D-g*C)*q;b[2]=(p*w-r*v+s*u)*q;b[3]=(-l*w+o*v-m*u)*q;b[4]=(-f*E+i*z-j*y)*q;b[5]=(c*E-e*z+g*y)*q;b[6]=(-n*w+r*t-s*B)*q;b[7]=(k*w-o*t+m*B)*q;b[8]=(f*D-h*z+j*x)*q;
b[9]=(-c*D+d*z-g*x)*q;b[10]=(n*v-p*t+s*A)*q;b[11]=(-k*v+l*t-m*A)*q;b[12]=(-f*C+h*y-i*x)*q;b[13]=(c*C-d*y+e*x)*q;b[14]=(-n*u+p*B-r*A)*q;b[15]=(k*u-l*B+o*A)*q;return b;
}

mat2.toRotationMat=function(a,b){b||(b=mat4.create());b[0]=a[0];b[1]=a[1];b[2]=a[2];b[3]=a[3];b[4]=a[4];b[5]=a[5];b[6]=a[6];b[7]=a[7];b[8]=a[8];b[9]=a[9];b[10]=a[10];b[11]=a[11];b[12]=0;b[13]=0;b[14]=0;b[15]=1;return b;
}

mat2.toInverseMat3=function(a,b){var c=a[0],d=a[1],e=a[2],g=a[4],f=a[5],h=a[6],i=a[8],j=a[9],k=a[10],l=k*f-h*j,o=-k*g+h*i,m=j*g-f*i,n=c*l+d*o+e*m;if(!n)return null;n=1/n;b||(b=mat3.create());b[0]=l*n;b[1]=(-k*d+e*j)*n;b[2]=(h*d-e*f)*n;b[3]=o*n;b[4]=(k*c-e*i)*n;b[5]=(-h*c+e*g)*n;b[6]=m*n;b[7]=(-j*c+d*i)*n;b[8]=(f*c-d*g)*n;return b;
}

/*
mat4.multiply = function(a, b, c) {
	c || (c = a);
	var d = a[0], e = a[1], g = a[2], f = a[3], h = a[4], i = a[5], j = a[6], k = a[7], l = a[8], o = a[9], m = a[10], n = a[11], p = a[12], r = a[13], s = a[14];
	a = a[15];
	var A = b[0], B = b[1], t = b[2], u = b[3], v = b[4], w = b[5], x = b[6], y = b[7], z = b[8], C = b[9], D = b[10], E = b[11], q = b[12], F = b[13], G = b[14];
	b = b[15];
	c[0] = A * d + B * h + t * l + u * p;
	c[1] = A * e + B * i + t * o + u * r;
	c[2] = A * g + B * j + t * m + u * s;
	c[3] = A * f + B * k + t * n + u * a;
	c[4] = v * d + w * h + x * l + y * p;
	c[5] = v * e + w * i + x * o + y * r;
	c[6] = v * g + w * j + x * m + y * s;
	c[7] = v * f + w * k + x * n + y * a;
	c[8] = z * d + C * h + D * l + E * p;
	c[9] = z * e + C * i + D * o + E * r;
	c[10] = z * g + C * j + D * m + E * s;
	c[11] = z * f + C * k + D * n + E * a;
	c[12] = q * d + F * h + G * l + b * p;
	c[13] = q * e + F * i + G * o + b * r;
	c[14] = q * g + F * j + G * m + b * s;
	c[15] = q * f + F * k + G * n + b * a;
	return c;
}
*/

mat2.multiply = function(a, b, c) {
	c || (c = a);
	var d = a[0], e = a[1], f = a[2];
	a = a[3];
	var g = b[0], h = b[1], i = b[2];
	b = b[3];
	c[0] = g * d + h * f;
	c[1] = g * e + h * a;
	c[2] = i * d + b * f;
	c[3] = i * e + b * a;
	return c;
}

//mat4.multiplyVec2 = function(a, b, c) {
//	c || (c = b);
//	var d = b[0], e = b[1];
//	b = b[2];
//	c[0] = a[0] * d + a[4] * e + a[8] * b + a[12];
//	c[1] = a[1] * d + a[5] * e + a[9] * b + a[13];
//	return c;
//}

mat2.multiplyVec2 = function(a, b, c) {
	c || (c = b);
	var d = b[0];
	b = b[1];
	c[0] = a[0] * d + a[2] * b;
	c[1] = a[1] * d + a[3] * b;
	return c;
}

/*
mat4.multiplyVec3 = function(a, b, c) {
	c || (c = b);
	var d = b[0], e = b[1];
	b = b[2];
	c[0] = a[0] * d + a[4] * e + a[8] * b + a[12];
	c[1] = a[1] * d + a[5] * e + a[9] * b + a[13];
	c[2] = a[2] * d + a[6] * e + a[10] * b + a[14];
	return c;
}
*/

/*
mat4.multiplyVec4 = function(a, b, c) {
	c || (c = b);
	var d = b[0], e = b[1], f = b[2];
	b = b[3];
	c[0] = a[0] * d + a[4] * e + a[8] * f + a[12] * b;
	c[1] = a[1] * d + a[5] * e + a[9] * f + a[13] * b;
	c[2] = a[2] * d + a[6] * e + a[10] * f + a[14] * b;
	c[3] = a[3] * d + a[7] * e + a[11] * f + a[15] * b;
	return c;
}
*/

//mat4
mat4.translate = function(a, b, c) {
	var d = b[0], e = b[1];
	b = b[2];
	if(!c || a == c) {
		a[12] = a[0] * d + a[4] * e + a[8] * b + a[12];
		a[13] = a[1] * d + a[5] * e + a[9] * b + a[13];
		a[14] = a[2] * d + a[6] * e + a[10] * b + a[14];
		a[15] = a[3] * d + a[7] * e + a[11] * b + a[15];
		return a;
	}
	var g = a[0], f = a[1], h = a[2], i = a[3], j = a[4], k = a[5], l = a[6], o = a[7], m = a[8], n = a[9], p = a[10], r = a[11];
	c[0] = g;
	c[1] = f;
	c[2] = h;
	c[3] = i;
	c[4] = j;
	c[5] = k;
	c[6] = l;
	c[7] = o;
	c[8] = m;
	c[9] = n;
	c[10] = p;
	c[11] = r;
	c[12] = g * d + j * e + m * b + a[12];
	c[13] = f * d + k * e + n * b + a[13];
	c[14] = h * d + l * e + p * b + a[14];
	c[15] = i * d + o * e + r * b + a[15];
	return c;
}

mat2.translate = function(a, b, c) {
	var d = b[0];
	b = b[1];
	if(!c || a == c) {
		a[2] = a[0] * d + a[2] * b;
		a[3] = a[1] * d + a[3] * b;
		return a;
	}
	var f = a[0], g = a[1];
	c[0] = f;
	c[1] = g;
	c[12] = g * d + j * e + m * b + a[12];
	c[13] = f * d + k * e + n * b + a[13];
	c[14] = h * d + l * e + p * b + a[14];
	c[15] = i * d + o * e + r * b + a[15];
	return c;
}

/*
mat4.scale = function(a, b, c) {
	var d = b[0], e = b[1];
	b = b[2];
	if (!c || a == c) {
		a[0] *= d;
		a[1] *= d;
		a[2] *= d;
		a[3] *= d;
		a[4] *= e;
		a[5] *= e;
		a[6] *= e;
		a[7] *= e;
		a[8] *= b;
		a[9] *= b;
		a[10] *= b;
		a[11] *= b;
		return a;
	}
	c[0] = a[0] * d;
	c[1] = a[1] * d;
	c[2] = a[2] * d;
	c[3] = a[3] * d;
	c[4] = a[4] * e;
	c[5] = a[5] * e;
	c[6] = a[6] * e;
	c[7] = a[7] * e;
	c[8] = a[8] * b;
	c[9] = a[9] * b;
	c[10] = a[10] * b;
	c[11] = a[11] * b;
	c[12] = a[12];
	c[13] = a[13];
	c[14] = a[14];
	c[15] = a[15];
	return c;
}
*/

mat2.scale = function(a, b, c) {
	var d = b[0];
	b = b[1];
	if (!c || a == c) {
		a[0] *= d;
		a[1] *= d;
		a[2] *= b;
		a[3] *= b;
		return a;
	}
	c[0] = a[0] * d;
	c[1] = a[1] * d;
	c[2] = a[2] * b;
	c[3] = a[3] * b;
	return c;
}

mat2.rotate=function(a,b,c,d){var e=c[0],g=c[1];c=c[2];var f=Math.sqrt(e*e+g*g+c*c);if(!f)return null;if(f!=1){f=1/f;e*=f;g*=f;c*=f}var h=Math.sin(b),i=Math.cos(b),j=1-i;b=a[0];f=a[1];var k=a[2],l=a[3],o=a[4],m=a[5],n=a[6],p=a[7],r=a[8],s=a[9],A=a[10],B=a[11],t=e*e*j+i,u=g*e*j+c*h,v=c*e*j-g*h,w=e*g*j-c*h,x=g*g*j+i,y=c*g*j+e*h,z=e*c*j+g*h;e=g*c*j-e*h;g=c*c*j+i;if(d){if(a!=d){d[12]=a[12];d[13]=a[13];d[14]=a[14];d[15]=a[15]}}else d=a;d[0]=b*t+o*u+r*v;d[1]=f*t+m*u+s*v;d[2]=k*t+n*u+A*v;d[3]=l*t+p*u+B*
v;d[4]=b*w+o*x+r*y;d[5]=f*w+m*x+s*y;d[6]=k*w+n*x+A*y;d[7]=l*w+p*x+B*y;d[8]=b*z+o*e+r*g;d[9]=f*z+m*e+s*g;d[10]=k*z+n*e+A*g;d[11]=l*z+p*e+B*g;return d;
}

mat2.rotateZ=function(a,b,c){var d=Math.sin(b);b=Math.cos(b);var e=a[0],g=a[1],f=a[2],h=a[3],i=a[4],j=a[5],k=a[6],l=a[7];if(c){if(a!=c){c[8]=a[8];c[9]=a[9];c[10]=a[10];c[11]=a[11];c[12]=a[12];c[13]=a[13];c[14]=a[14];c[15]=a[15]}}else c=a;c[0]=e*b+i*d;c[1]=g*b+j*d;c[2]=f*b+k*d;c[3]=h*b+l*d;c[4]=e*-d+i*b;c[5]=g*-d+j*b;c[6]=f*-d+k*b;c[7]=h*-d+l*b;return c;
}

mat2.frustum=function(a,b,c,d,e,g,f){f||(f=mat4.create());var h=b-a,i=d-c,j=g-e;f[0]=e*2/h;f[1]=0;f[2]=0;f[3]=0;f[4]=0;f[5]=e*2/i;f[6]=0;f[7]=0;f[8]=(b+a)/h;f[9]=(d+c)/i;f[10]=-(g+e)/j;f[11]=-1;f[12]=0;f[13]=0;f[14]=-(g*e*2)/j;f[15]=0;return f;
}

mat2.perspective=function(a,b,c,d,e){a=c*Math.tan(a*Math.PI/360);b=a*b;return mat4.frustum(-b,b,-a,a,c,d,e);
}

mat2.ortho=function(a,b,c,d,e,g,f){f||(f=mat4.create());var h=b-a,i=d-c,j=g-e;f[0]=2/h;f[1]=0;f[2]=0;f[3]=0;f[4]=0;f[5]=2/i;f[6]=0;f[7]=0;f[8]=0;f[9]=0;f[10]=-2/j;f[11]=0;f[12]=-(a+b)/h;f[13]=-(d+c)/i;f[14]=-(g+e)/j;f[15]=1;return f;
}

mat4.lookAt=function(a,b,c,d){d||(d=mat4.create());var e=a[0],g=a[1];a=a[2];var f=c[0],h=c[1],i=c[2];c=b[1];var j=b[2];if(e==b[0]&&g==c&&a==j)return mat4.identity(d);var k,l,o,m;c=e-b[0];j=g-b[1];b=a-b[2];m=1/Math.sqrt(c*c+j*j+b*b);c*=m;j*=m;b*=m;k=h*b-i*j;i=i*c-f*b;f=f*j-h*c;if(m=Math.sqrt(k*k+i*i+f*f)){m=1/m;k*=m;i*=m;f*=m}else f=i=k=0;h=j*f-b*i;l=b*k-c*f;o=c*i-j*k;if(m=Math.sqrt(h*h+l*l+o*o)){m=1/m;h*=m;l*=m;o*=m}else o=l=h=0;d[0]=k;d[1]=h;d[2]=c;d[3]=0;d[4]=i;d[5]=l;d[6]=j;d[7]=0;d[8]=f;d[9]=
o;d[10]=b;d[11]=0;d[12]=-(k*e+i*g+f*a);d[13]=-(h*e+l*g+o*a);d[14]=-(c*e+j*g+b*a);d[15]=1;return d;
}