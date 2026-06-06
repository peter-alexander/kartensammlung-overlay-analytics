const B=`https://${process.env.EASYNAME_HTTP_HOST||'tiles.radlobby.at'}/Gehsteigbreiten/Testdaten`;
const P='work/sidewalk-widths/';
const W=[[0,'< 1.2 m',null,1.2],[1,'1.2 m - 1.5 m',1.2,1.5],[2,'1.5 m - 2.0 m',1.5,2],[3,'2.0 m - 3.0 m',2,3],[4,'> 3.0 m',3,null]].map(([id,label,min,max])=>({id,label,min,max}));
export default {method:'fmzk_guided_sis_multiband_cross_v1',crs:{metric:'EPSG:31256',pmtiles:'EPSG:4326'},paths:{workDir:P