# ملفات COMSOL - التوثيق | COMSOL Files Documentation

## نظرة عامة
يحتوي هذا المجلد على نماذج COMSOL Multiphysics للمحاكاة العددية (CFD/FSI) 
للتحقق من صحة النموذج التحليلي.

## ملفات المعلمات

### Parameters.txt (الإعداد الافتراضي)
```
H=14cm, L=16cm, D=12cm, A=30°, f=1Hz, Φ=40°, N=3 rays
```

### Parameters_last.txt (الإعداد النهائي المُعتمد)
```
L=14.4cm, D=9cm, hmin=3.4cm, hmax=17.4cm
A=30°, f=1Hz, Φ=50°, N=2 rays, ncycle=6
```

### Fin_parameters.txt (زعنفة كبيرة الحجم)
```
H=17.5cm, L=20cm, D=14.6cm, A=15°, f=1Hz, N=3 rays
سماكة الرابط=3mm, العرض=0.6
```

## ملفات النماذج (.mph)

### النماذج الأساسية
| الملف | الوصف | الحالة |
|-------|-------|--------|
| `unduFin-FSI.mph` | نموذج FSI الأساسي (3 أشعة) | الإصدار الأولي |
| `unduFin-FSI_v1.mph` إلى `v4.mph` | تطويرات تدريجية | مراحل تطوير |
| `unduFin-FSI_statinary.mph` | تحليل استاتيكي | مكتمل |
| `unduFin.mph` | نموذج هيكلي فقط بدون تدفق | مكتمل |

### نماذج case studies
| الملف | الوصف |
|-------|-------|
| `unduFin-case1.mph` | الحالة الأولى |
| `unduFin-case2.mph` | الحالة الثانية |

### نماذج أخرى
| الملف | الوصف |
|-------|-------|
| `moving_mesh_laminar_wave.mph` | اختبار توليد الأمواج بشبكة متحركة |
| `mechanism_submerged_in_fluid.mph` | آلية مغمورة في مائع |
| `project1.mph`, `project2.mph` | نماذج مبكرة |
| `Fin_mechansim.mph` | نموذج آلية الزعنفة |

## مجلد Fin/ (نماذج فرعية)

### الإصدارات الأخيرة (last versions/)
⭐ **هذه هي النماذج النهائية المعتمدة**
- `2rayunduFin-FSI_Cd.mph` — نموذج 2 شعاع لحساب معامل السحب
- `3rayunduFin-FSI_Cd.mph` — نموذج 3 أشعة لحساب معامل السحب
- `Cd results.txt` — بيانات معامل السحب (3005 خطوة زمنية)

### حالات الدراسة
- `2Ray/Case1_2Ray_phi0/` — شعاعان بدون إزاحة طورية (Φ=0°)
- `2Ray/Case1_2Ray_phi50/` — شعاعان بإزاحة طورية (Φ=50°)

## بيانات النتائج (Cd results.txt)
```
الأعمدة: Time [s] | Sway Force [N] | Surge Force [N] | Heave Force [N] | Cd [-]
الخطوات الزمنية: 3005
المدة: 0 → 3.05 ثانية
```

## ملاحظات للاستخدام المستقبلي
1. النماذج تحتاج COMSOL 5.x أو أحدث
2. تأكد من استيراد المعلمات من ملفات `.txt` المناسبة
3. النماذج في `last versions/` هي الأكثر نضجاً للمقارنة مع النتائج التحليلية
4. ملف `Cd results.txt` يمكن قراءته مباشرة في MATLAB بـ `readmatrix()`
