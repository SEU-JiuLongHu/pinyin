#include "HMM.h"
#include "PinyinConverter.h"
//#define _check_vector_size(x) if (x.size()>x.capacity()*0.5) x.reserve(2*x.size());
HMM::HMM()
{
	
	init_pinyin2chars_table();
}

HMM::HMM(string corpus)
{
	init_pinyin2chars_table();
	loadCorpus(corpus);
}
HMM::~HMM()
{
}

bool HMM::loadCorpus(string corpus)
{
	ifstream fin(corpus);
	if (fin)
	{
		string wordUtf8,pinyin;
		wstring wordUnicode;
		uint32_t cnt;
		while (fin >> wordUtf8)
		{
			wordUnicode = PinyinConverter::Utf8ToUnicode(wordUtf8);
			vector<string> ss(wordUnicode.size());
			for (size_t i = 0; i < wordUnicode.size(); i++)
			{
				fin >> pinyin;
				ss[i] = pinyin;
			}
			fin >> cnt;
			for (size_t i = 0; i < wordUnicode.size(); i++)
			{
				add_char(wordUnicode[i], ss[i],cnt);
				if (i)
					add_chars(wordUnicode[i - 1], ss[i - 1], wordUnicode[i], ss[i], cnt);
			}
		}
		return true;
	}
	else
		return false;
}
id HMM::get_charid(wchar_t ch)
{
	auto iter= char2id.find(ch);
	if (iter == char2id.end())
	{
		id charid = char2id.size();
		char2id[ch] = charid;
		id2char[charid] = ch;
		return charid;
	}
	else
		return (*iter).second;
}
id HMM::get_pyid(string py)
{
	auto iter = py2id.find(py);
	if (iter == py2id.end())
	{
		id pyid = py2id.size();
		py2id[py] = pyid;
		id2py[pyid] = py;
		return pyid;
	}
	else
		return (*iter).second;
}
//½¨Á¢Æ´ÒôidÓë×ÖidÖ®¼äµÄÓ³Éä¹ØÏµ
void HMM::insert_pychar_relations(string pinyin, wstring chars)
{
	id pyid=get_pyid(pinyin);
	if (pyid2charidlst.find(pyid) == pyid2charidlst.end())
		pyid2charidlst[pyid] = id_list();
	for (auto &ch : chars )
	{
		id charid = get_charid(ch);
		pyid2charidlst[pyid].push_back(charid);
	}

}

void HMM::add_char(wchar_t character,string pinyin, type::count cnt)
{
	id ch = get_charid(character);
	id py = get_pyid(pinyin);
	char_freq[Key(py,ch)] += cnt;
	py_freq[py] += cnt;

	emit_out[ch] += cnt;
	emit_to[make_pair(ch, py)] += cnt;
}

void HMM::add_chars(wchar_t character1,string pinyin1,wchar_t character2,string pinyin2,type::count cnt)
{
	id ch1 = get_charid(character1),ch2=get_charid(character2);
	id py1 = get_pyid(pinyin1), py2 = get_pyid(pinyin2);
	Key key1(py1, ch1), key2(py2, ch2);
	transf_to[make_pair(key1, key2)] += cnt;
	transf_out[key1] += cnt;
}

vector<Res> HMM::query(vector<string> pinyins,uint32_t topk)
{
	matrix.clear();
	matrix.resize(pinyins.size());
	observation curobser, preobser;
	state prestat, curstat;
	probability prob;
	for (size_t i=0;i<pinyins.size();i++)
	{
		curobser = py2id[pinyins[i]];
		auto charidlst = pyid2charidlst[curobser];
		size_t charsize = charidlst.size();
		matrix[i].resize(charsize);

		for (size_t j = 0; j < charsize; j++)
		{
			curstat = charidlst[j];
#ifdef DEBUG
			matrix[i][j] = MatrixNode(id2char[curstat], curstat);
			Key curkey(curobser,id2py[curobser], curstat,id2char[curstat]);
#else // DEBUG
			matrix[i][j] = MatrixNode(curstat);
			Key curkey(curobser, curstat);
#endif			
			if (!i)
			{
				matrix[i][j].logp = log(1.0*(char_freq[curkey]+1) / (py_freq[curobser]+charsize));
			}
			else
			{
				for (size_t k = 0; k < matrix[i - 1].size(); k++)
				{
					prestat = matrix[i - 1][k].stat;
#ifdef DEBUG
					Key prekey(preobser, id2py[preobser], prestat, id2char[prestat]);
#else // DEBUG
					Key prekey(preobser, prestat);
#endif							
					/*int test0 = transf_to[make_pair(prekey, curkey)];
					int test1 = transf_out[prekey];
					probability p1 = matrix[i - 1][k].logp,p2= log(1.0*(transf_to[make_pair(prekey, curkey)] + 1) / (transf_out[prekey] + charsize));*/

					prob = matrix[i - 1][k].logp + log(1.0*(transf_to[make_pair(prekey, curkey)] + 1) / (transf_out[prekey] + charsize))+log(1.0*(emit_to[make_pair(curstat, curobser)]+1) / (emit_out[curstat]+1));
					if (prob > matrix[i][j].logp)
					{
						matrix[i][j].logp = prob;
						matrix[i][j].pre = k;
					}
				}
			}
		}
		preobser = curobser;
	}
	return solve(topk);
}
vector<Res> HMM::solve(uint32_t topk)
{
	static auto cmp = [](const MatrixNode &a, const MatrixNode &b) {
		return a.logp > b.logp;
	};
	assert(matrix.size()>0);
	uint32_t L = matrix.size() - 1;
	sort(matrix[L].begin(), matrix[L].end(), cmp);
	vector<Res> records;
	for (int i = 0; i < min(topk,matrix[L].size()); i++)
	{
		int cur=i;
		Res record(matrix[L][i].logp);
		for (int j = L; j >=0 ; j--)
		{
			record.str.push_back(id2char[matrix[j][cur].stat]);
			cur = matrix[j][cur].pre;
		}
		reverse(record.str.begin(), record.str.end());
		records.push_back(record);
	}
	return records;
}
void HMM::init_pinyin2chars_table()
{
	for (auto iter = pinyin2chars.begin(); iter != pinyin2chars.end(); iter++)
		insert_pychar_relations((*iter).first,(*iter).second);	
}
const map<string, wstring> HMM::pinyin2chars = {
	{ "a",L"àÄ°¡åHï¹ß¹ºÇ°¢ëç" },
	{ "ai",L"°¨°¥V‹Ü°§è¨š±êqËB‘°ñLðgšGƒù°®ïÍÌ‡BàÉ°£”±ˆì•lŸs°¦­a°¬ìa´o°¤°ŠèPŠâœâæXæÈƒvÙŒ°}ŠÖ‡†éu³}òIö°‰¹Üt”²êÓ°«êi°ª°©œÜëB´Ì@…¥†ì\àÈ‰a²}Þß×cá{íÁ³v•á°¯øµKžG°­‘¹Û×r…Ù½i" },
	{ "an",L"°²±Vˆ¥°°Øtƒ‡ëˆ«q^†±ŒåðÆÚÏÁOÞîóì”ÉŽë@áí›¡õc°·…{°³÷öí™±Q°¸ÈCˆù“°¶ˆÝ°µä@åBéœñKîO´UÇIùg°±ÑsÄWÕY‹jV†HÈ€‹Fí°´âÖÖOèñ††Ès¹ãï§••³§¯uØÛû" },
	{ "ang",L"Œì°¹–‹óa°ºáZál…n°»" },
	{ "ao",L"âÚ÷¡à»…–À†õ“³ˆŠWŽS°ÀÊTª‡‡Ìéá°Ä—`ø€á®C­H°Ã‹®ŠSÂK“ýåÛúqåòˆòü’UæÁü‘RÛêæñ‰¥÷é nŸÑ’jÎ‚öËëJ°¼Ö’ö—EÏùÂOéO‹‹ÆbñúÞÖ°½°Â‡ÆÒ\°Á°¿°¾Ö“ÝE´x" },
	{ "ba",L"°ÏŠ‚ŠB…©«X°Å÷Éþx°Ôˆ¢Ú•°Ñ‰Î°Ë³Fy°ÓášþwRÁjâZü––Â™ñ°Æ°Î°Ð¼“ÝRØ^Ôy†^°È°Éå±–[Ò†õNÃ_”îÙ÷ˆöÑÁT°Êá±ôÎQ°ÌžßïTõE÷„°Ö°Ç°ÍÜØ°jÝÃ‰°Õá—ˆzÍM ã°Ò" },
	{ "bai",L"°Þ”[°Û–à°×ËbêþªW»Ÿ°Ý²®ßÂ°Ø°Ú¸q»“Þã’…Ù”Òo†h½]”¡ív°Ü’“°Ù" },
	{ "ban",L"Š”–®°çÚæŒêÃR¶t°ß°äÞkE–DÛà½O°áîC»O°â•LâkÞnÎŒ°ëé›ñ­îÓÛAñ£ô‘°ìœ°ã[ì‡»{”‘ÎZˆÐ°éˆm°è°àÒƒÑ—°åøX°ã°í”Ê°ê°æ­šô²" },
	{ "bang",L"°îÖr½‰¿R°ð°õ°ôžI°ñŽÍ°ó³‰«g—”Í{°÷ÏÅŽ°ÝòÍK‰Y’ÊäºíD¼’²°öŽÀóo°øˆÈ¶œ°ùÎM°ï“s°òß™ˆ  ¥æ^" },
	{ "bao",L"ˆóìsŒp±¡Ýá¾‹Ç˜ˆç±©æßÊ}ÆØÙèÍdèt„ô«’ÅÚ•Þƒ˜–¢å²±ªõÀ°ü„ƒÑfAÒJï–õUöµødóbÞ‹~¸±¬ÅÙÆÙìÒð±é–°ú±¢ñhŒšÜ°û±¨±£±«è˜ìdŒ—·‘óŽý_ï’±§°ý±¦±¤øRÌ™ÐˆãEŒ‡°þÙ…ñÙ±¥Ë" },
	{ "bei",L"¬i±°¹u÷¹–È±ºðÇ±®ˆ¢ì‹ØÃ° †h±ÛèE±³±¸þw±²ÕRãm¬DªNíx—À—f‚ËÊ‚pÒoíÕ•K‚³±¶¶F ÍØ‘v±·àf—“ã£ä^Æp±¯Ý…ËöÍ±±±´“d—G±¹±­¼L±»ÝíßÂ±µùlÚýƒFÚé ´‚äÈiñØócÍ“–{ÝK" },
	{ "ben",L"êÚŠMÏn›y“àªŠœ`†Ï—LÝ™ Äº»ÙSï¼èM±½±¼Ÿø‚–ÛÎ±¿’ÙžÇåQ±¾—ñÛÐßGÁ" },
	{ "beng",L"ˆ©±ÁÐ±ÄÈ’²°öÛM¬e‚õ¬aìž¿‡éa±À±Å¯n±ÂÈEê´éGˆÈ¾XåAçaŠR‰lµpàÔ½l±Ã“gßJ" },
	{ "bi",L"–©†ô°nú‡±Ë±Í–Äæ¾±ÛÜK·Kó±ÒˆfáùÉœæÔÀVé[®n—éÚFíSâØ±Ê‰ý°zøpÚPª‹·÷ÃÚŸ•¯w®wÎ“ÁXØ„µ–Ý©ÜL»zß›±ÏÃZ±ÔÓvÃØ»õÏ¯HèE§¯RŒù›ñÔ‚¿ð{œ í{Á‘¿oääîé±Ýé]ös±È–ŠÔvå¨ÞµÛ‹Ø±Ü¼ž±ÕõIßÙèµ±Ì”Àó÷ôxí@Ã^ŽÆésÌYØP†žŸÎêÚî¯ŽÅôÅ±ÚãG±×ü„àˆê\Äb±Øò±ÓÙÂÙSì‹ÒK±Ð÷Â±ÇW±É«ÓÍšÏô°šÈñƒŠ`õmåöóÙóëñEÆƒØ°ªŒP…ñÐ‹ˆãàŠÜêš¾aº`¹v”è±Ù×’®…ŠŒY¸“Ògã¹Û~±ÖåþÈ]ïõŒÂœü±ÑùS˜[Â›a±ÆïßÁÝÉæqç@ûGÙC¹P÷”é\±Î‰ªúz—a–aŠË" },
	{ "bian",L"öýÒŒ ¤±æ•còùÅX±è‰äËx×ƒñÛçÂÓS¹Þq±á·HâíÞlÛÍöc±éÌ¯Vª ìÔ¼D±ãÆß…›M±ßãêì™ÙHöbÈqß„±äÅŒáŠóÖøuíÜÜÐÞg¾œ±Þ±ç«fí¾»e“Oñ¹æQ±åØÒªpîYÞp±â¾Ž’\érØPß›±à" },
	{ "biao",L"™~˜Ë÷BÃ ƒšïlÒFªY÷§±ìždñÑ‹›ès·…œýïkïRïÚïjì®ñ¦ûïðïn÷ÔÕ•æô±ë´‚œWÊEì© gæ»‚l‰wŸÏì­±íålÄróTìá±êè¼Ž¼ÅAÙ™ï[Ö€“¿òŠË‘çSóQ" },
	{ "bie",L"õ¿±ð÷M±ïü‚„eÆƒ°T±îÍr•…ñÇaÖÏh±ñ–Äý–ÒXÌ‹" },
	{ "bin",L"¬ž”PÚSó‰žIß“ôWžMÙÏÀ_÷Þš›î ±÷óë÷ÌžÙšÓŸØh±õÙfžlÄœ±öšñ±ò÷ÆáÙ±ó­pìEÙeïÙçÍ—ÃéÄƒ†™‰óxè\±ôéëçã" },
	{ "bing",L"±ú±þ‚§—Š–â—€±û¸pã°SÆÁêv±ý²¡·’±øðV˜ðÚìh‚ìÍsšê¬V’m•½Ž±ù°R•m•\KŒ}Úûïž±üÕ@éÄ–Þ±}Ù÷™‰âT²¢ŽÕí@ŒÎ’ò½lãu·AÞðÆuŽð‚víSõm" },
	{ "bo",L"ñgópŠ‚ÑJñ•±¡ò’°²²÷Q²µöÒq²·ë¢²«“ÜðG²³°hãK´BØÃÈ• éñCº~³jÒT þÊXœÀÆ…²°íçàR­“ƒ`ÍoõÀ™q–Â²®ŠÈ`„ƒã\Ü@óŽ“²§¼\™ØéÞÃJ ¦ñFÔy²­›ÂõÛ²±œ_­”Ÿ¹ÆÇµR Ý²¤îàêþÑBœ”‡¥ð¾ÙññA°l—K†\²ªâÄõN°Øäc‰®ðo÷ˆÐ“è}’©¼žÃ`Øm‘Åô¤à£‡hÊN¶z²¯â“ùPÒU×LõE²¬‚Nªt²¨´‘ƒkænÆt²¦²¥²©ómÀõËÅ‡ŒX°ã²£°þéD¹ÌY±C²´" },
	{ "bu",L"ÆÒm²¿²¼’p‹åÍŽï²·²ºÞKÑõ³„Ïß²øGÛYê³ðX¹rÝ•–¿…ùÕcÇ[E’ÃêÎ²¸²½šiŽ~îÐ²¹ûQ…Ä±¤îßÑaùLŒ “äà^º^ÉžðJ²À²¾ªŽ¶âb²¶ê†â˜Šçˆ¶šh²»" },
	{ "ca",L"²Áàêµg²ð”cßníå‡Í" },
	{ "cai",L"š¶²Ã²Ä†’’ñ¾ZØ”²Æ²È—²Ç²Ê²ÉŠéÛPˆÆ‚š¿n²ÂŒu²ËÀu²Å²Ì‚Æ“H" },
	{ "can",L"‚ð…¢è²²ÑÐQ²Ò‹ìïŠ²Ï‡A |þp²Î‘KÎ]²ôœ’ï{šˆ÷“ß÷õ²ÍÐT‹ÛÖƒ…üo•ü Nåî‘”öŸ‘L“·‘M“½ôÓ†Ð²ÐæîöY…¤ÓËL²Óò‰ÛŠ" },
	{ "cang",L"‚á²Õ…MÙ‰Ø÷û]²Ô²×Å“ª‚}²ØÊiº[ƒûÈúIœæ²ÖžPÉnè†Ï@¨" },
	{ "cao",L"²Ü‘¨ç[æàÐ²Ûýó©äîòxó“Ùåø²Ý²Ù•ùÆHÜ³ÒG‚ó²Úô½‘FÉ˜" },
	{ "ce",L"¸žÈYÇRŽú¹kƒÔ¹‹²à‚È…‹®‚Å²â’‘²áÉƒ²Þýv²ß”˜œy‰x»‘Š¹ZâüÈmºu" },
	{ "cen",L"ä¹—q³•…¢²Îßá¯¸’" },
	{ "ceng",L"²ãò¸}àá•û™IòšÔø²äŒÓ" },
	{ "cha",L"¿’·Ñ–àêŠg²îñÃèd²ê²çÇNÆOãâ²ï²ææ\ã˜“ ïïâªé¶æ±²í’Ké«ìx¼pïÊšÛ‚÷ðlÔûâÇÔˆ“c„x–Ë²é²ëè¾ÃPéßåš²ì²åâO’Qˆ“†â’¼‚²ÓÉ²²èÅ‘Œð—^ÅaÔŒ¶g" },
	{ "chai",L"âOðûÐƒýb²î®› åÏŠÙ­Óµ}‡Ð†¶ò²ƒŠ’K²ðÆ²ñîÎ²ò" },
	{ "chan",L"‡Áã@ÀW„•Ægä²ûêèÀsÏMÙææöð’²öÝƒdéˆª†äý„­‰Ê‰²ø²üçP„}ÊrŽf›ºå¤Õ¸žeºoÒ—ìøš´ÖÝÛ³ƒ—{õðƒ{ÏsÃˆáp•Càš˜^²ô²ù×¶U¨ó¸êUÒR‡Ïâã„iPƒ]I™ÙŸžåîµ¥ÒcÆB®b†Î²÷”â´v×€åñ²óæ¿éK®aàž“·äaÐŸ“½í]žŽÊœµ‡c×‹äiÚÆ‘ÔâÜ†®ŽÂ‘Ï×ã²úž¬“˜Þ{‹ÈîèÕS²õ¿Cƒ§ŸíÀAÕ~Òb—ÑgÀpïâ”v" },
	{ "chang",L"ÌÈ²ýŸ…•³æÏÄcáäâê³¨³¡Ÿ³£ü÷löðÕkéL«`ãÑ•˜ØöêOS®^äœCå_÷•ã®ƒ¸ÈOêÆ¬díoÄq Ã›éMæ½ä–ÝÅ®D²þ³«³©¬„…”‚tÛËƒ”ˆöÉÑë©Ñm‡L³¤³¦®˜Ï^öK¬ çL³¢³¥è ³§—–³ª‰jƒYé‹‡ŸÜÉ" },
	{ "chao",L"´Â³¬Ÿ·ÞCÁV³±„à¿žÚˆìÌÔN˜ÈÖšŸq³®³¯ûž¾K³°³­€RÚ}»}à}êË£³²â÷¾bŽlŽz“¼ Ÿ³´±|³³Óe„¤ñé¸Jžânü{½Ëü…™ù" },
	{ "che",L"¼‚e³·…ãŸLØÛå³ºÖ‚®Šb†q …³¸åø“ÝŸEÂs„ïÍ’îJíºÀÞŠ³¶ÇpŸc³²uÔaÜ‡“F³¹³Œ³ß³µ" },
	{ "chen",L"·Q¯M³¿³•—²Ø÷“ZúmýZ¾DÇ_ƒ¡Ÿ‹ÚßÖnå·³Ã³¾ÇkÒr³Æ¯’Ü•êÏIÖR‰mé´’×ÆëËl•æáp³Ó”³Â³À–×ÚfýYÞÓí×Ú’Ù•‚áÔH³½êJÞ×¯„¿bû‰‰}ÚÈÊc‰ö˜¹Õ€´~²_—FàÁ‡¸×â\è¡³ÁÙo³ÄÛ{™Â³¼åŒö³³»" },
	{ "cheng",L"·QÚX™fˆá³ÏòGõ““£Ç^³Ò‚ èKîª¸V›„“Î’Þ‘ÍÀ˜çd‚Dîñ³ÅòÉ³³ÆŒk¿BŠ¿‰S¬A³Ð¬šÛô«žèßÚW˜üØ©±›Æ–b³Çšé¯·œ‘r³Ó—¢¹f›Õõ¨ ªîõ³È³Ë¬b˜û½†®—ç•ð‰³Í’¬îd‚tŽñ‘³Î³Ñ™rëó‘~žj»ñÎ›“çpÕ\wÃ”˜Œ³Ìä…Ï|àá—¼·Ãw æjêpžsèÇÊ¢œË³Éìl—–³ÊÍBàJòr\êÉ" },
	{ "chi",L"–okˆk¶ß³ÞÞ‹øTÛy¹}ÚdšlüJ¼YÙP³ÚßgžÃ‘dýXð·«”~Á‹Ö–šIÇK¯øï†Õvß³Íh›n³Õó×ñÝ»ŒñYÜ¯‘y…µ›‚ÍN„„Ã’òˆ‰‘´Óó¤’xŠw³Ø°VõØ³×ÙÑßoúuø|ã‰ñ¡…q¸‡ü[ÑDŸUÛLâÁ’»Ú—ßê³ÖóøÄS…äÃqÔ ÃL³áà´áÜÖs‡iãM“¤Úm³â…Õš“ÃnŸë³àÑEò¿¯b³Üë†³ÙàÍ½‚ÔWÚpß†ùAI…h¹xßt‚sÂ]Ú†³ß²lÜÝ¹M³Û†Ëå~Û­˜»’³ãÂ@³Ý®Eí÷ßWœ‰p¯€íôãrôù ô‘JÛæ„È¯vÂB“¹÷Îù–«„ÐÞŒýcë·ÝBuÎyù`ÕBæÊÑlµo³Ô’Lð„" },
	{ "chong",L"³å«–Û ï¥ô©¯\‘oã¿“_³ç³è¾…Üûƒ†ü×ÛŒÓ¿ÐnÏx³äô¾›_ÖÖÎuã|ê™ ‚“›ÁZâçrÖØŒ™³æµr‚ò†Á˜¶›ÒˆÃÁˆÑ~" },
	{ "chou",L"þr–„±T×‰¼—õ–ƒ½[³òáh²ƒ ¶® áOƒ‰õ\ÛS³ð â‹á³êšŽÜPÅW×žÔ—‡œ³óºN³ï³ñßc×p—¹“oã°Ñn…Á³í»IÙ±‘ÀËg®‡³îñ¬³ë×‡³ô‚¸‚G¾I–ä³é³ìŽÎ™„±y°{àüJáböÅ ß bÇ“Öaá~ô{ël" },
	{ "chu",L"Øa÷íÂ^ñÒ Ë”™”ßÓcÆcÂaèÆãI½IË ´¦ÔxÕ‘³÷™»éÓ|ºX‚m’}³û´¥éË´¡˜ýƒÛ»‚âŽÐÐE‡bžã´¢³úÊxÐóˆÇ“ªØŒ‹ƒÉeýs™úƒ¦¬GÌŽ¬`ÉZÜX„IØ¡Éâð³ùƒšb­lèúàsœä¸aÛU‘Ãõé¸e´¤äzërúR—ÆØX³üòÜÄ•ýi³þãÀ“¹ç©™sŒçézN³õ³ýÚn³øÛH³öµA´£‘Aúž" },
	{ "chuai",L"õßþŽšIà¨ÄuàÜëúÞõ´§þŒ" },
	{ "chuan",L"å×Åxë°â¶´©„”ÝŽîËº@´®Ûw•ÄçÝâAÄxúE´­ô­šö´¬´«‚÷´ªªk¬ÇFÙi´¨šN‰@«[ƒb‡ù" },
	{ "chuang",L"êJÊ[¯w·™„k‡l´¯´²‡è ¡´}„V‚üâë‚}´°´´ §í´±„y“œ–S„“„€ —¸R™H´³Ä€¨" },
	{ "chui",L"´¸é³Ç”ôD´·Úï“€ÄD–û´¶¹Š‚…æm×µé¢îqˆ§´qý—åN´µ´¹‡ù" },
	{ "chun",L"õžÙƒÆXÝ»ù‡É”™šác´½œ·Ã–~ÈN•«ê_•IÉO›Ì´¿Èo˜‡¼ƒùœÃ‹´º´¼‹a‚¤Ûw¹—²Qþ”¬tÝòí˜êœ÷ðÈåTöj‰@ Æ´¾´Àƒb´»" },
	{ "chuo",L"´Â‹S¾Yà¨ïß´Á¿žýwæ—áQåYŠÆýpöºÚ}šf›íŠÅÞu“óèqßOär´‡¾b·‹Cê¡‡Ç×ºÄJÝzõÖ" },
	{ "ci",L"ôÒËj×ÈŽú²îÞoóq¿WŠœð@–Ü´É‚½°rŽãÞi´Æ´ÃÜëËF´Ä´Ë–æÃh«y«R´ÅèÍyâ‘ˆôôÙ–cÚežBÆ˜–Ÿ«uÜùú\´Í´ÊµQÚ…‹½a‹ãÆú]øy´ÎÎˆÕ›õJ„pý€ÈW…è´ÈßÚðËÞeËÅìôÇ„´Ç®N–²´ÌÕ†Ô~ódÏ…ï“ˆˆÙn" },
	{ "cong",L"ÄÖÆ‰Ê[´Ó˜BÏZ™ß|^œYÙz˜º³Ÿ‘m˜ÚäÈÂ‡—Œ´Ð ´ÑÙ{ËqžšÀò^çW¾‘…²²jæõ Qƒ´ÔŸtŒQ‡èÂ”ÜÊÉÂŒèÈò•¾çýè®Šæ˜âÀSS•›´Òß´ÏŸÐºbò‹­B‘FÕp¾" },
	{ "cou",L"ê£œ´ÕÝëíé¨" },
	{ "cu",L"´Ù×äâ§¿qÚ‚¯|Úuû›Óc´×î•ÃÛq¯•õ¡´ÖŠÆÛUû‚û€‡m…aÜAÊIÕKŠÅ‹{üy´ØÛn‘–ÛcÚ…éãÊPÝý’Ûû„õ¾åeáÞõí" },
	{ "cuan",L"´ÛŽmÇˆÜf™«ìà´ÚägŸäß¥·‰”€”eïéÒ{´Üš–è‰¸ZÔÜÙàºx”x™ç" },
	{ "cui",L"õ˜§Ãœª‰»‚ßýÄƒ‘N´áŒ¾\ºÀŠÄ‹´àîx¿\û¼´Ý´äÚ~ƒþè­¸Wé¯QyŸÕ¿…Ãy´âö¿´Þ°„ÒPéÁçJ—½†Ÿë¥‰…´…´ßã²¬X´ãÝÍ‚yË¥ŸnÄ›ÁŒ" },
	{ "cun",L"´ç´å¶×›–üââñåÄ~…¼’Žß—¸€„Y´æ´»vÛZ‰–" },
	{ "cuo",L"ðû´ëØÈóq¿W„v´ìðî´íái´èÉxÌ‘ÉcÊPßuõãåe„z´éà™äSÒPï±Üg±‘¬›áAÇs‰èý€ïóÌßHáÏàŸûzõºÇuÕ´ê— ëâ" },
	{ "da",L"“ÒÛQÔz´ñßQœÍ÷°ÇQèNÞÇÜJæpñ×óÎ¸—í³…A‰¡Þ‡æ§ËR‘„®}šùæ]‡ÀJý‘ËþÏƒ÷²í^žØˆ™ý“´ó´ï®†´ðÇE‡}Þ‰âòðãÞ…ß_´ò„‘´î…ößÕÁe…ì±oàª" },
	{ "dai",L"‘·ÍfÛº‰´õÜ–ølþåÊÞaÝDÜÚ±¬xß°K´ûá·Ž‘ÔrÏEšù÷ì´÷Å•µ¡´ý´ú´ø´ö–±¿DÛFünì^ÙJ´ô´óªyìOñ~žŽÜ¤ß¾Î}çª´þçé•Î´ùŽ¡ææŽ§…¦½H´üÊO¹yˆ‚Òy" },
	{ "dan",L"ŽÀWíñáG‚„Ú›“úà¢…S‘žµ­ñdðZÙÙƒd†›ð…ÕQÊ¯¯Dµ¦á]Ä‘µ®å£êæµ¨µ«‡dÒ—‹[µ£îFº„ÖñšÂnµ¯µ©ó‡ìKŠl›X…gÂ›Øéü^“ÛÚàµ¤‡nµ§ÒRµ¢‡·Ül’b´‡~·žÏ€ÄEÉ…µ¥×¼“ÚÍž†Îø}µ¬ñõürš——‘„·ééóìÙyµ° ý­®X„éÓ”†²Ùœ°QÝÌàülÓg³Nðã„[Ÿí°D–½ð÷µªÒb¶VÐy«m" },
	{ "dang",L"ÊŽ«šÉ‹ C™nº‚¹YÇµ±×•ñÉ“õ‹PÚèKÕ‚«ˆWšë±Uµ²ÚÔ®”ºš²^µ³µDÝÐµµêWµ´´XÌoüh™éßT¤­c×[­Tå´ë‹ƒ}žªÅ™Ï}‡ŽÛÊÒdîõ”†‰³®Gˆ›í¸" },
	{ "dao",L"Ås·RëZµ”µ¿µ½’Òµ·Ðp÷ô€µÀâáß¶µºÈKŒ§‰»¶\¹|ƒ‰ c‡‹ÅŽW±IëŽÎìÂRµ¾ÁŸÍ@ë®ê‰µ¸àüµ¼Ü„ôîu–]á’­”Fµ»Ñn½rìâµ¹“v™|ëIµ¶µÁÐm" },
	{ "de",L"‡NµØzÃŽµÄï½Ôµ×µÃ›ú¿œåuµÂ—‘" },
	{ "dei",L"‡NµÃ" },
	{ "deng",L"ç‹ŸôàâïëëQØOµÊ¯µÅíãê­Å˜­O³ÎµÈáØ™žƒ\µÆ¸~µÇô£µË‰œà‡‹¿µÉÓR" },
	{ "di",L"’†ŽRˆkÊH¶EÛyÆlÚd¼e™˜Ûqãd‘dµ×÷¾µîE‹X‰y˜µµÓf~—\¯FíÆ‚±—bˆhˆ¹Þž›‚ƒCÚ®wæ·µÑµÚµÝÚÐÙáé¦àÖ®Sˆ¯ƒ™ÓhµÔµÒóž’FµÙËy’ãÆmâKô†«Z•AœvôÆµÌ”³íÚ†vÔgµÖÇ… ¹íLÄVêÌáñV±ƒµÛ“W…}›ÁïáÛ—–š‰„‡”˜NµØÕœµÍÏEßróƒµÎµÄ†¬ûMª–mØµ†OŸb‰—…àÝ¶tÃJµÕó†Kå~ØpÖBÐ”ÞË‹Û¡êsµÜµÎ[µÏÇœßfêëºa«Ÿ“ŸÓ]Â‚ŠD´”ÛæÛ‡èÜíû‚dÝBÊLçCœìµÐÊO´YÚh¾†ì{µÞ" },
	{ "dia",L"àÇ" },
	{ "dian",L"Ês‚Ù‰«ëŠ˜•¯tµáµê‡Ãâšµëµç˜ˆÛãîµâµà¯’†ˆ‹L—Ï”¥ŽoÛ†µäîFñ²ŠHØ¼”„Ôa°dò›‘úçè˜ëñ°…ŽÑµíîäµæ‰|ŽpÕ¬Uµå´ô¡É_áÛµîý‚µãµèîŒÄHücõÚµßµé™AµìÚçŠû”“ÍŸ" },
	{ "diao",L"µïÉ‰âyäHµñõ ³HõM¸uöôûbµõëµòèS„aïM’FõõµöŒÅµ÷ÓŽËyÍ@ÍqÄñ·–‡ÝUô†ï¢ÕA¬h¸Lîö®œ@µóåcøJÚwÕ{šôã“ŠPù@½rážÀ’ü—‹àä”²føB¬tùmµð¯šµô" },
	{ "die",L"LÛ’¡Þég•iõÞšÛ‘ä«µøÕ™íCç“•èU˜GÔeÏH¯Aµüõ]µûëºñóà©†AµþÃ]Â½xÛìð¬˜›À„ÞµúöøÒBšŠµùéPÜ¦ÂWÚg±yÅŽÎH š¯BÛLµýÛ@õÚölÆ|›u†—Å\Ñ±‚–»" },
	{ "ding",L"ì´Oï}¶¤åVàËYí”ç–î®ÆJ¶©Šc—ÅÈbð—Øêîú’ðîr³G¶¦ëë¶£çà¶¥íÖâÎ_á”à¤´ŽŠäbÓ†ü‡ñôôúðÛ}ìw¶§ÍB¶¡¶¢¶¨" },
	{ "diu",L"îûäAïMG¶ª" },
	{ "dong",L"ëØá¼Õ‰âº—òLßËù…ð´ŠàÊÄL¶­¶«žúªJ–¶³¶¬ˆÄüŠ’œ‚”Ô˜–|Þ“žöCá´ÛíÐhñŽë±ÎXdƒP¶®šæ“_íÏ‹Ùœ§¹š‘ãƒö›ò¶²¶¯Ç‡•kŠŸ½p‰’„ç¶±¶´ëËõ[ÁÆ{„Ó¶°ëš¸•" },
	{ "dou",L"ò½ž^¶¼¶ºñ¼Ã–óû„r‚J¶·¶»…Êô`ÝúðLšÃ¶Áî×ô^×xƒÃêL¸]àKôZ¶¸›Ãêh¶¶äWÇWð„E–’†táH¶µ²fôYôa¶¹ékâ^”Ô—u™X" },
	{ "du",L"÷Ç¶¼…XÚG­{™³˜ÌÇTÒl¶ÃØKÐCÜ¶‹óá`¸]¶½äÂ… íbÒe÷òêA•’ªšÙ€Š¶À¬o°ÓG¶ÉÎ–ê^¶Êíž^òy¶¾í~èoŽª”¾ê•¶Ä Ù„†æN¶ÙÑtåƒ¶Â×˜×x¶Áì|èü¶Çë¹¶Èš˜„‹óÆåLà½šœ¶Å¶ŠîD•¤Î}ºV„E ©ó¼ütÔŒ¶¿¶ÆÕà" },
	{ "duan",L"ÂZ¶ËÜYÑƒ¿EÄa‚Ç‹e‰Fº@óýìÑ¶Ì„Œ¶Ð”à¶Í¶ÎÈ˜š¬å‘´V¾„»f¶ÏæH¬‡é²" },
	{ "dui",L"×mŒ¦¯yêŒ¶ØíÔ‘‡í­ø‹Ëc‘»¶ÒçŽ“€ƒ¶îX–€žSˆŒæmç…‰[¶Ô´qžwê ×B¶ÓŠZåTÅž}í¡žA¶Ñïæµq½˜ƒµˆÍ" },
	{ "dun",L"‰Ý¶×¶Ø˜ú÷ìÀõ»ÜO“æßqãçÞš¶Ý‰•ÜHç…¶Úí»ŽÝâg“ÇÄ]ˆdïæ˜Jíïíâ‘‡ª–¯Ûv¶ÙÎP¶Õ—çŽŸõ Ô‡¶ÛîDªÄR¶Ö´]ò—¶Ü" },
	{ "duo",L"’–’—Ü€¶ß¶ä–šÛT¶á–ÃôDñÖ”£”¦ˆ‘ßÍ—ÙõâêyˆÊõy¶âñWÛGï˜„m†Æ‰š…¼ç¶–m‡šßá‹sÆ–ÚrèÞŒ¹–ú„„Øy–ª¶å¾E‰™îìŠb¯k¶ç–\¶à‰ï”ŸÍÔêwð™“¶ÈèIãûü„‹”­‹µÜo¶èÛF›k¶Þñjàùz¶æ‘†ŠZ¶ã¶éšÇÔq‡¾„AëD" },
	{ "e",L"¶ó×†´dêqÔ›ÓFÓžùZ“t‡fŠŠˆº—¿™ÄÍL†s³j°x…vˆì³bù[ßÀ‡êÎY†¡Þˆ³SºÛëÚÌŽþ¶íð_Sî€ÑÆŠã”Aî~¶ô«Å¶Ý­åíîPòFãµæ¹ÜÃîOèyÊ‚³r†H…ÅœŠëñ¶ïÖ@àöˆñÛÑâ…¶õÜ—¶òãÕ…Ù†@‚­ù˜¶ö³X¬Ø`›áöù“épôŠÝàˆ×™šdýL¶ìò¦ïÉŠŽðIï°âeŠ´Šâô‰¶ð¶ñé‘¬c‘öØ¬åŠ×FðÊ”šß]÷{†‘ý|ý…êiƒišxŒï°¢ötŒß±“ß{kø‘–•ÝQ“~ÕMÌFéî¶ë„þ¶êä~¶î" },
	{ "en",L"ŠCÝì¶÷ÞôWàÅßíŸ¸" },
	{ "er",L"¶ûñ“Ù¦Ý[êzçíÐ^ÂxƒºŒª‹èÙEÙ@ð¹åÇËnp˜ÞîïõbßƒÂYê—–éÃsŒ©Ú·¡¶ü„n –¶øX¶þ¶úÑLÇHó“›˜¶ùó’öÜ–êðDø¶ýš¾ãsƒ¹†„ÞWëX…þÔ õ" },
	{ "fa",L"ÙHÛÒ–ì·¦°l·§‘“ÜšøóŒ¸ŸÁPáNÆž·£óŠíÀ˜ì¬m·©²XÊ†Š‘žž·¨›oáw·¤·¢åz¯VÁU‚ëéyËt·¥Uá" },
	{ "fan",L"”ó·µïcÅxœt˜õÅwž’»o®‰ïxËXïˆ—÷ÓŒé·¬÷Y·­ŠïâCÒT”õ·¶‹ÑÝGšï‰“ í·°Þ¬çx·¹…K¾u¢î²·¸–¯¹ „å·³Øœ·´·®¹DŠi·«·ºï‰·ªë¶ÅtõìÞÀ’Bšø¿œ·±·²‡h»Oáë„G··—¡‹Ë–iú‹á¦Ÿ©Ä‡·¯µ\Ï›±F‹Ì¹B„FJÜèóÁ€‘ŒÞN­[ìÜž~" },
	{ "fang",L"¼ °ØÎ·Å·Âô³·ÄáÝœE·ÀîÕ·Áøh·¼ó„‚·½­œ·ÃÚ“±}ˆÚÏˆªÔLèÊ±fÚú•XöÐ”ëúJ·¾µp·¿ô™•PÍKâ[·»" },
	{ "fei",L"ì]·ÇªU—’·ÈÏn¹Aìqö­–FÊˆ•›öEáôöîŒÐÕu·ÐŠóì³½EåúÊ†ü–·ÏïyÈQ·ÊUÑpü”ç³·Ëã­êŠ·Î¯XÎN…ŠòW„|•Õð[¾pìéÑqŠôŠO·ÍÙM·ÌëèïÐçšÒU·Ñó‘é¼ñI˜ìôä·É•hÊ„äÇ–ÉÜÀžOóõ‚n°CÃ^ïwòãòa·ÆÃd™¶‰–{œdðò" },
	{ "fen",L"ÄÞMëƒ²b·ÞèûñOÉkÁi·Ü—±÷÷ŽËŠ^ÃRŠ}Êˆö÷–DîC¼SÍ`ŸþŸøƒf·ÔÁ‰ð·Ó·Ûk·ÚôšžÇÐv·Õ·ß³WØk•SŽŒ·à—rÓŸˆe¼Šä—ƒÀ¸jÈ†–Œ·Ý·Ø‘Í_–BåªŠ°ä·ÙØrðiå¯èMüv·×·Ö·ÒÁ‚÷aˆbâpñBëV™JüRü‹¶løXÙÇŒð‰žçã" },
	{ "feng",L"Ÿ‘‚ª·åœtçQÙˆ·â„O¥í¿müKìbºA´^‰â·èøP·æ—÷·é“žØNoŒ›œ˜¿p½ §·ãÙºãã±`›ÍÅ‚ƒt·ï È·ìÝ×ÃTÚRñTˆ©ÛºÖS·ë›høiŠ~™l¯‚Ÿuähœ½’¸ŸÔßô—QÌtà•ØS¬SÅ}ž–·çïp·áæ‘ïLÇl·ä·ížÐªhˆùÒƒ·î„KÐIøL„N·ê®g" },
	{ "fo",L"¦–ˆu·ð" },
	{ "fou",L"Švžäë€À¼€ÀŒø]ó¾Æ]·ñ" },
	{ "fu",L"Ùë®i¸¢ÜÞ½•ß«s¯ž±GÍb–ó¼”žÞ¹AãVTØ“¸¡öûŠï·òöv}áUæâÍ|Ì’Ð•øDŠ•Æ…–ŽU‹cïOäæ—ÓÒi·ûÝ³ÈƒùfØfŽ“›š¸¬¸®—­‘Ê¹[¼J‚¾¸½“áå˜º…Ðuî\õÆâaÉ’¸£·öÖÍkõvíhÒLáK¿`ÁÝÊ¸¦¼›¸¤¸¯¸·øWÑ‡ûŸ¸±«cÇXèõ¹r·÷ˆ}¸©åõïTÖD¸µNãR®w›^ívˆ¡Û®íê¸«·óºŠÂ—Úç¨Œ á¥®tÈiòó’ÑÓ‡‹DøqõÃ¸²¸¶õH¶@Æ]ôfß»ìðˆóÛ~ŸJæÚ‡`Í½n¸º†b·ðµyƒå„_î·ˆŽÝ•¸À®}Žˆ¶OÎl»™¸­ù›õV½E·ü‚¿¸¿ß‘–¢ê‚Ùx¸´’½¸¸¸ªüFÍ—·ùïûíÉ·Jáœôï·ýÞÔíëÜ½–´òÝ¸¥·þ·öò¶Ò„í‚ÙŽŠmÇCç¦üA¸³·ôÜò¸§àMð¥âö¬MÑ}Ð“ÝPÔcÙM‚a¸à~¸°Ý—·øó‘¸¼·ú¸¾‚YøIŸr­oÄw˜_–Áå‡êç”êÝoÜÀÅ€Ê…ò³Q¸¨¾”›Lñ€¸¹ƒì¸cÙìòð·õ›Š¹…¸»Ãi½šÁJT" },
	{ "ga",L"¸ìáÞÎîÅ¼Ð¸Áåm‡QÜˆÙ¤æØ¿§ßÈŠAê¸¸ÂæÙôpÔþ«V" },
	{ "gai",L"ÇD„÷½w[êà˜£¸Ç¸Ä®„Ø¤Ô“ã¸Æâ}˜¢“©ëB_­yØdŠ¡¸Åµ‹ÙWYÙ^•|¸Ã¸ÈÛòÉw½æÈ‘Úëê®–qæYì„øà@½i" },
	{ "gan",L"«\Ús”Ôl±Y«qðáãïÛá¸É¸Ñ¸Ó°‘ÜÕ˜o¸ÎÞ|‚‰q»ˆäÆó_Î’I¶’öxŒÀx¸Ê¸ÌŒ¼åœØJºTí·º•ôûÇ¬ç¤ž¸ÚC¸Ï÷h™gÍH—UŒ¿›NâFÞÏ¸Òôvß¦¹m¸Ðl÷ ƒ÷¸ËêºÆQÐr½CéÏ„QøN¹CÚMŽÖä÷¸Í" },
	{ "gang",L"¿¸¾V¸Û’âŸ€„‚í°¸Üƒé³Mî@¸Øêl¸×‘ßœÏ ±¸ÙæsâG‘Þ¸ÔŸƒÀ Œùñþþb¸Úóà˜þhÀ“ˆÕ Â¸Õî¸´Lˆþä“¯I—ž è¸Ö" },
	{ "gao",L"°wçÉÕa¸Þ¸ßÅV…Ì¸äï¯ÁoË›z¸á»ízÚ¾ªˆ•±éÂê½µ†úkðp‰ù¿cú„Æ¸â˜°™Ró{¶ŽéÀ™²üŽ²G·X˜‚µ‡¸æÞ»¶JÛ¬ä†Ç¸àæ€Øº™¸ã¸å¸Ý¹lœõ" },
	{ "ge",L"ªn¸ìÒÙ¸óøœèÞ‘áÝ‘¸õ¿©ÛÁ¸÷…ÃíkéxñËíu¸è÷…’MíR‘ëÛÙ¸ç¸ëæü ·”R’š çíÑô´¸êÖgàÃœð–qëõ¸ñÅZ¸ôÓkéwõsôŸïÓ˜†¸íòZ¶…¸ï÷Àë¡…Ï¸éÖY”šÍÞPÑ\¸ÇæŠÔ†ÃIÆŒùBæk¹wømØîÄ—ò¢¸öÍxøwºÏ ³„ý‚€™ ¸ðò´Éw¸î†þì‘îMàõi†ñö¸òª˜ÜªØªãt¼v" },
	{ "gei",L"¸ø½o" },
	{ "gen",L"Ý¢¸ù“^Ø¨¸ú“jßçôÞƒ" },
	{ "geng",L"¹¢ùˆŸ‰’ù¸ü¹¡ÙsÐÏûfç®aîi¹£¾±½Ž¾•œÈ@âÙ½b’ªßì¿K®uàQ¸þˆíƒói—ÔÇcõ†¸ý½c„jÁ}öá¸û›Êy" },
	{ "gong",L"…CÁ‡ëÅ¹§†ßÚCýŠ–íô„ÜpÞÃ¹±¹°‘E’¹«…š¼k¹­ŸËò¼‰b¹¨ƒÅäUãŒmýœ|´bì–Ý\ö¡¹©¹¬Óyºìó•¹®¼t¹ª¹¥–r¹²çî…@¹¤ºT¹¯ÚMØ•¹¦†yŽ³" },
	{ "gou",L"âh¹¶“Â¹´ÂV¹¸ƒÚ˜‹ÑåÜ¾—œÏ¹³íxØþÂU÷¸á¸Ý@Ôì°¹º«vóÑ’]¹µˆxæÅÚ¸ã^Æa¹·Ô_Š¥‰òÆ™ÓMº¹¹ÂTëgèÛ‚×Ím¹»ÐêíŸµ¾äçÃØxóôÙ“k" },
	{ "gu",L"êöˆØ¶™±W¹ËîÜÝLÉu†g†åßÉ°›—›åd¹Çâ’’MÙZ¹¾‰à¹Áãé˜€†f÷½õY¼ÖúX·Y¹Æ˜bù]âÝž†˜ðÀÔb¹‡‚ïì±Ë[ÝM¹Í¹¼™O¹Éôþü‰¹½¼ÒÁlöñòÁ’_³‘žkî¹¹Êëû¹ÌáÄƒóøŸ‚¹Âéï†Ø¹Ã¹¿ðó¸š¹À¹ÈÆ‚¹Ä°–Ú¬ØÅ›}Hêô¿S›üµÁBÐMëºHèôõýð³ƒlöAÍvð –¾ÝÔžJí‹²Œ½î™¹Åðkî­ïÀ" },
	{ "gua",L"„œòmÚ´†§Ÿ…ÔŸÛ|¹Ó¹ÒèéˆqßÉÕ ã”Æ‚šO…³ëÒ¹ÎøŽŸ°ð»¹Ñè’ìäTÁGÀ¨ƒÖïW½\ÁL¾ Úo„Ž¹ÐØÔ¹Ï" },
	{ "guai",L"¹Ô‰øÁLÞâ‡ˆ–Ês“–¡¹Õ¹Ö¹y…¨" },
	{ "guan",L"¹× ƒäÊð^²›Þè¹ÜÉF¹á¹Ù÷b¹Úµ¹ÞØž°H÷¤æšêPÓ^îÂOŠþ›Œ‘×¬gÝ¸¹àè…¯p›ý÷}ßk‘TûX­À•…j“¥ÜIév¹`µeÒ‹¯¾]¡­e¹ßêKå]ðÙ²•ÓQ ëqÂÚÅo¹Û¹Ø˜ÀÙÄ¸AÝ„‚š¯¹ÝøAÈXöŠ" },
	{ "guang",L"ƒZë×¹â‚UïŸD›²ž»ˆ’•VèæžÖÅQÆšÞ‚«EÚ‡ßÛ¹äžÕ“ÑÝ_™¤áîžÓ¹ãã «‡Š­üU³qÅS" },
	{ "gui",L"—ËŽQµƒ•Q”¯ÒŽÆ—Ž`²Z”‹ˆ’˜³Š¹¹ê‹‚Ú‘é|¹îæ£ôkŽ@‹¥æ²n˜­ê{öÙ„¥¹ôóþÄ„ÀL¿þ…QŽë¹ðÙF°I«•²zÉ}ð§÷¬ø_@…‘¹è˜²÷Z¶W¹éÒ^íW“±ªgœÄêÐÎš¹K˜œšð“Êßž³u™Æ¹çøWÍŠ¹ñ­YêÁõqØÛ‹¾ØÐÚbàF„£Ãv‘Ü‰—Î™™÷iœˆ–_šwèí¹æø`¹íÈ²¸Y¹óôhŽ¢­„òo…TâÑ¹ë¹ò¹åëvÏj¹ìå³ºlÑO™uiý”w·˜ÔŽÓm¹ï”Š™Í" },
	{ "gun",L"çµÊFÐ–Ö²O¨Ùòõ…—œ¾iõP­eåKÝ¹ö¬g±šL¹õ¾É€íÞöçÑr¹÷" },
	{ "guo",L"œàþßÃÇ‘Âƒ†J¯†ë½ÙåœuÝ{”šÄBÄNÛö«†©¹ú‘I¾[‡ó²Ž¹üèJÊbòä‡î¼@â£é¤ã‡ëðR˜¡Ž½‡Hß^å¹øòå¹ù‰˜ñøÄsÎÐ—ëÏXðŸˆåÞâX¹ýâu‡øºl›ý“‡ñáÆ¹ûÑx" },
	{ "ha",L"ãx¸òîþÎrŠoÏº¹þŠU" },
	{ "hai",L"áVëÜõ°½wñ”º¡ß€¿Èï™º¤º¥ðŽàË…õº§º¦ŸQ‡¯»¹º¢º£ñ›ì†ãà@ºÙ" },
	{ "han",L"ØEº­÷ýÊG•Âˆ¥Í”²ën”º·ê\°y¬HŒåÎLú[Î‘˜oÞþäwô_…î¹bº³Î‡öîu’I—êR•~…{ìÊí™ñU›þÎKŠÎº«º²º¨åØJœõAŸß—cº¹òAºªº®ºµ†c„T†iÇt‚þäI›ÛÍHòÀ›¿ädša—UîhÖ››N‹©º¶ínâFþ[º°º±é\º¸ªRº¯êÏÝÕ®]›ÈÚõ‡•ÃQhºº¶îMþ\ãÛ—ßÈÌkñHº¬ñüøNå«º´ò¥º©•ˆ±Ž" },
	{ "hang",L"î@º¼ÍaãìôŒº½”ã¹VñþÏïÐÐç¬Þ†ß’Ø˜º»Šs½WÆf¸‘¿Ô" },
	{ "hao",L"•a¸äª‚àãhå©Â|ò«ºÃºÑ»ƒŸž®•±ºÄ•µˆÝï¸hÞ¶ºÁö‚ê»°€ºÅÏ–œB•¼ð©°‚¶mË^ª|‡_•Ø†Sòº•‰àÆºÂ…ëœéº¿ºÆ°…°ˆ»DËAî—×qÂG—·Ì—Ì–Ër‚Û‡sæ€Æ’ºÀšµÕ’àzå°º¾" },
	{ "he",L"’uàÀÔXÞºÕºKºÉðgŸŒÔZéxûi´ºÑ—æ÷…ÛÀºÎ…ôŸZýLÏÅ†ÛæüÎ˜’š ç±AíHÝ °œ¸ºÇ†YºÒûSúQºÊëaŒyÚ­´EàA²ˆ†ÛÖôŸý˜PîÁšBŠéu¶…ºÖ _–­Ÿ¿¼† ôçý†ÏšÏ½œzºØ°FùŸÔ†ý[úKÙRêÂèY”—ºÐãFãØÔò¢ÂGºÏÐŽºÈ»t‡˜ºÌªCºÔ eÖyº×±BêHºËîMšÎûÒ‡ºÍÈMØ€ºÓðšŸÀï¼v" },
	{ "hei",L"ü\àË¦ºÚºÙ" },
	{ "hen",L"ä‡ŒºÝÔ‹ºÛì•ºÜºÞ’‹" },
	{ "heng",L"ä™Mºàºãžî†‘›êˆýaçñÃ†èUûaèìÃtÐÐºâºßø’Š¬ºáÞ¿" },
	{ "hong",L"éb‡«³…Þ®é•ºêºä†ßºæãüåÈˆ“EØFœ‚Å|…šüZìºç…Æ¸së”ãµ³{ŽcºåŒâ½“È‡ºì›ÍŸp…ÔÝ¦ØAëŸ¼tž¿ø™ØD¹ÙêÙä»ŽÁ‡ô\“Ðºéšºëô„Œf›ÂoŠ¼ÁŠ›K›Äˆ˜~ãpÞ°•{Ý“¸fÞZäU…·œ|Á«YÓÚ§«aé{ãÈÆyäfÀ€é—¼‡ºèŠk¼˜âv†yºC" },
	{ "hou",L"ºòÔýJõ`àCö×È‰åËðf…Ë÷c³@àjÄDÂJã²TðúØ_›•ºðÂFÜ©…éºíóóááŽ«÷¿ºï êºñºóˆ‹÷æAô×ö\ºî" },
	{ "hu",L"óË—üºn¹”‘ò‰Öø‡úCò®¹}í’•Uì²Ì•ìï•÷ºý‰Ø…IðÉ›R‘ñá²ù]“‡ÊSÔS›~ÙüðÀº÷ÆS½œœû®@œWÏ·ŽðbÐk‹|ußü×oºüGàñí_†Ø‘õúK÷sžCºû‡Pã±»¦ºô½`”NÌ‘ïºø›´âïìÎºùž€ºÍ‡©»£»¤ìæçú˜«õúà‚µCÆ~ä°Êd»§îgm­•ìÃ¿eôE»¥ëŒŸW÷½›Z—ý‘ôúX»¢[ºþ…OƒêÓ{·‚ºõÄŠëiÖ—õ­äï‚såt‹­åŽÄéõÐí’_“ªö{ÈL•Oºú÷ŸŒŒØmˆ~šXŠýù–œX‹¬¿SHÎ™Y›üûIŠ¯»¡ŸÚºöºË–ô–øUð­éÎÜ â©ìè†¼ð×æL‡F" },
	{ "hua",L"´hüX™Š»¨»ªÝ{ÀEÌsÌf»®ÈA†®“Õ ­Lò‘Ô’±Ói»©“®‹½ÉJ¼@»¬‹Ã‰þ…ÅèëÕ–åkâEän“çÖœ×fæè»í–Î”Š£‹O‡WÅp˜¥çf“Š»­Êyîüô‰»°®‹ú†Æ_„»¯ªœ»«˜å" },
	{ "huai",L"»±Ñ‘Ìx»´žxÑœÌ|‘Ñ»µ‘¯»³‰Ä™Æ»®ÝÂjõ×†F»²" },
	{ "huan",L"ÉVÃK»»çÙóO÷ß”kšZ‡ÈÛ¼†¾ŠJ»¼ùJïÌ»À´öéØo½b¯ˆÝÈœo í¬~âµÞSš÷äñ‘×»¹»¾šgà øb»½±šëf²o¼]‘¤Œ~Žw×Ÿ¨éIûXØhË»¸ÀQÈP“Q`­hûqêX˜¬»¿×’»¶à÷»ÃÝkÑÛ¨ß€­’±Ø}»·å¾†¿²`öd»ÂžðÇBêaõŒä¡åÕ‡õæD½Œ»Á»ºˆâªBÈëq™ö¾ØŽß§öZä½ŒAÁv‚—hèG" },
	{ "huang",L"U»ÊéB›²äÒ‰E»ÌëÁÐYÅŠñ¥è«ðc–M°ŒÔ…Öe»ÐåØ›Róò˜R»Äéí‹»ËÚ‡™¤‹h»ÉŸºòböü»Çå–»ÆŽxœê•Íª‰Ÿ÷UŒr†ÅƒÆ»ÏáåömæwäêüS‚µ“NÖWÈˆðçuó¨ŠN·k»Å»Ñžê˜nŸìúŠ¬‰•sÚò»Í»È°¢»Î" },
	{ "hui",L"Ñ‹çÀŸ˜•þËD»â‹^ä«×eÌsÂP›‘»Ù›i‡éš§ßÜó³¬qÍ »áÚ»²»Ò›x»ã‡j³|ØYÊ]’’‡Gã„»Ýš«ÒË™ŒáÙVö™…R­_Ôœ¡‡v»å»×íW“Öž`‰™÷âŸ@¯»ÛßÔçž²Nª›Ð„ÁšÀDŒ@—ÛçiýIèíÂEh×M‡‚»ÓêTj»æ˜ž»Ü»ØÞ’Ý™uåçËCýH•Á»ÔëD»ß‡¤Üîî_ÍY»ÕÍzÝxÆU‹ÏìuÖM»ÚŸÌlƒaà¹Ž¹ˆHïH“]ÀLÞ¥îœ×w—ò»à S•Ÿ Z{o™BÚ¶ŸF@ê_À£ŸCÔÕdò³Ò^ž¾êÍ™bÎšƒª»Ö¢…¡Œ“ãÄí£Üö™mœóßD‡ß»Þ‘Îçõ™®í}­g»ä·xÁ™œ“þiŒà¢õt×fº_‘}š•Šîµ˜‰Ä¯`ðd²~¶éŒä§" },
	{ "hun",L"—p—yãÔ»ìœ†—•»ë•eý@ÞF²Jš‰ÓoŸ[ða±d²EÊMùœ¡ÕŸ»èçõ¬q»ê¿Œ¾iäã»çî‚»é›÷ÈÇ‡õðQ‚“’ä“]é’ˆâ‚[ùŸkÚ»‹GâÆ¸" },
	{ "huo",L"ïÁÞ½¶„œ­ÄNâ·»õì[èZ»ò»ðØ›»ôtòdŸZ´žµœ¶åxÉ^²‘ïì·‚Öf—ë„Š…ü±nôr™Š•ëŠ_ëoÂhó¶‚ižmÅŸ‡—ß˜àëžCñëâ€«@”N”üß«…¿ìáîØ»îØå›[»óè‡É»í’»ÄsŠ£’îå»ö°\ºÍÚoéXÈuÅG»ï²ˆ»ñ" },
	{ "ji",L"ë}½o»÷³€·bÜÁ°nÆïÁaÞU”úèWún„©é®†ÀÚ|¼§Ÿdçáî¿Ûˆ†æ°U¼ªUÓfò±óÅÆnïWÅU³ÛÔ˜Û¶Sˆô»ûÎa¼¾ì¥¼­‘¼É°uÛEîÜQ´‰çˆÆÚ…¯ÜùÌIýRï|¼°í‡ãšŠjMêª»ý×IìV¼ÅöSc¼´¼±‰IùH­DÑ_È—“ôñ¤•¸·}Ï…Ó‹úWŠ¸Ý‹št“Ä¼¸¿ƒ“VßóŽóêéç¼ÈöÝÂf»ù”ûôn¯séêÄlœ–ÆëÊDÁYØÀÙŠÕ‚ëYÒQê÷˜O™›ƒÎ³ß´½Y¼¢‚ÂÅ ä¼½½åð¢ØG¶P¼²—mÏµëu¼ÊØ½·eÙ}™Cú„Þ‘¢^ìPÕHÓs…uŒïçg‘¼ÌRÊmê«œg–ˆÈWúaÑw”ªÞá›ù¼Ì¼·çÜóÇÆˆ­^ÜuûAŒ¨óKÛ”Ò‰²¹UàBÁb¼oÙ¥¹œˆ…Ëj¼¨Ûp‡÷q¼¡÷‚·]‰€¼¬»üÒHð‡˜Šì´õÕ—ùºsÏlýW¼¼á§Üe¼Ææ›‹¼µê÷ä™vûn™‹Ûeèi™Ã¼Á¼ËÜ¸ì“Ìzä©¼¿ÙÊøK„W˜œ¼¥ÆI¼£ÎŽ¼‰ƒ_õŸô‚”Œ¼_æ÷õJõÒ²]»þš©Òˆ™o­u¼Ã·m¼©¾@Ž×öaˆj¼¯å‰ÂcØÞºuÀMÌn»úÞªØ¢ßâýVöêíZÓ“Ï„ë|Ì~T¼®åæô‡÷D¼ÂåìÕ‘¼ºýTåZž†‘ÕÀ^¼ÀêåÓJ¼¦…hëH¸ø„Z¼³ØCòT¼¹¶I‡\ŽNÝð“Ø¼Ç‰JÆ–Æä·I¸ï¼Æ¼¤ïú÷C”Dã‚áÕÃôßö«’í¶®‚Ž©ªE¯Ó]•ÌŠ Úl¼Í»øËE„ˆßÒ¿]÷Ù¼ÄÓ›™W¼»Âö›Ûa›D¼«ù¼¶¾ƒ×^Ÿ" },
	{ "jia",L"ãxÛOŒ_“ü¹kÝç’z¼ÞàP’~ïØÙ¤Ž·˜–æ‰¼ÙÙZîa’¶áµê©˜k¼ÏÄ`›Ñîòû“¼Ö’S¼ÐõÊñÊ¼Ü”Ï…­ƒr¼ÚðèòÌ—kóÕÃˆ]î]¼ÑÇv˜\ëÎ—ÝË¼Òñ{«w¼Ûðø”ðýš¹ªmÂæˆ®Ø†›vÂ_ÇÑŠAÑW¼Õâ›åÈøm”ÐØÅ¼Oäe¼Ôñ˜Û£ä¤Í™xí¢‹TÐ®ùG¼×‚í¼Î ÇãeØjôÂ†kçìðšªo‘æ¼Ø†í¼Óï¼Ý" },
	{ "jian",L"û{éfçZª\˜cübâ]ƒc‚¡ídÛ`²væZ¼ùå€ê§ïµâJ¼üæI½¢ÆDÚÉ¿ ¼þ¼÷²€žR×vÅ[¼øËuöx¼î„‡¾àîÞY×t—ä›–‚ßÀo‚kÝÑ½€Ï•Š§”W™zê¯Ú{™Ò›×žˆåÀ‘âìyš–çõÂ¼ó“B”ƒ€™‘ŸÒûy¼âµM’³ã‹ÒŠ¼ä‰AäÕä’„§¼öégƒïÊ—Ÿæ„¦ÅžÝóÈ‚øZæG±O½§åbå¿ˆÔÔdå¹Ùvöž„ªšž¬‚í[ÒZÕ‘ìÈ“‚›”ð¼åçÌuóÈ²Rû|ûx•©ÈG“b´š¬{á_èÅÑI«…öäâV¼í¼ã½¥¼ôåsº†¼ñ½¨ÙÔ¿VôC”s½£¼àÊz÷œ÷µèƒ¼éœ—„Éë¦ÏM¾}Ú™òq—gúYžh‰¤Ø]½´Š¦´DÓS¼õ½¤êùôå˜öñJ¹aç™  íKêðœ\¼òé¥èaö¼ìè—ëU‡Ø½¡Ë]´–’þç‰™ZÙ`Ò}¼ýŽ¥õÝ»W¼û¼ð¼G½¦ÒMüû…ž“ìëìÈ…ðeíú²{Ì‚è~ŽÔœpùp¹{—ÊÞöÖGŒ{¼ç¼áÇ³¼ïÊ`º]ØbÂžÓVÒOÚÙèB˜ÙžŒÖˆñÐä[è{ÀO¼ß—ß«lðÏ¼ú„«Àwæ~ör¼êèb¼ë¼æðT¼è" },
	{ "jiang",L"®{…G˜ªí\÷FÇ¿ä®áuŠ\‰á@îŽÑH½³‰Öv½µãÜ™Äv®–ºç™ºçÖ½¬ÏQán{Œ¢¼Tôø“°Ê@½¯÷šñð½©ÁžôÝ½ªÜü½²š™ÊYŠÚêñ®Ÿ–tÎ…‰¬ÀP“Àç­ËKŠX½®„ß½{½«x½°Ö˜íä½±½´‚×½­™^Èwª„" },
	{ "jiao",L"½Ìõo”©‘xž«­dÞBjú„õäÐ‹É·p³C½ÂÚŠÜú†ý½¶þ™]ôé‰½½“è”‡ÅT•Ý¸‹†û½ÏÙÕÄ‰°½ÅÛ]ÆLÀU½Ðƒ‚½ÆÀq”Ò¼iófÚˆÒ™ð¨kŸ”œ©½g‹ùƒeÙ]ÜF•w”œ÷R½¸½»á†ïœ’¹‡EàÝÄ_ùaª—Ð£‹´‡Uáè«„ž•Ÿ÷½¼½¹½ÉŽBÝ^½ÄÙ®Æ›ÓXöÞžìˆ’›ç€°‰½º‹Ð×_Ïf²½À¹R HËŠ½Ê“¼œò½È½ÇãqõÓ½ËÓŠ™ËÞIðÔþõ´½¿½Á„à·•‡„á½ºŠÒŸ}½Î•¯‘¢úŒòÔ”º½Í“×æù×K½ÑÄzÌ—òœ½ÃÞØÜ´…ÓÏtæ¯½¾½·„¤”¼ë¸¾õ" },
	{ "jie",L"™w½Þ³V½ã‡‚ÍÈÔ‘½ÚÓnŒô„o—ô½ë½ÝÎf½×˜‹ˆêÐVŸ®—¸¼½Ó½Õ½æµ@ÀTËÇÕm‹‘™ÃÚµì“½ÒÙÊÎa…mÞ×æOðÜò»ôÉ”T’ù ÏÑ›”â½Ù¬p‚Ü›­‹d°Xô‚½â˜HÏÕ]ŽY†‡õ^ìŒŒî¹œïà®½Ø’÷ÖŠæÝ•MÂcÛOÃöÚÐ|ÑKÍßâ½ßÓ“š²æ¼“ø…Ã½Û½éÐwŽ^½äíÙ˜mëe÷º…‚Œˆû¯^Þ—½Ü½ìôÏ˜ŒÃã]Ú¦Í„½ê«d½àù™½Y“ƒ½ÖœœÇ}„Â”O½å‹}½Ô½e¼Ò„f¼Ûð·MŽÑ½çÉ•Ñ\ºÍŽªEëAò¡½áÛdŽà„g˜PÚl½è®vîRïN—A‚í›Œ¨àµƒËÐèî¯CŠo†Ö‹m" },
	{ "jin",L"½ñ˜cû²›½þŒƒ„BûvÇM½óš½ü¬’äWêá½ùÚB†‚‚BÇžáŽƒ½÷ßM”Ü„ÅâY«ƒ³\ŠúÝÀ½ýçÆƒq•x½ø¼Žƒ»‡ž­\µÓb¾¡µ‰ƒH¾oéÈË|âËÓPèªÖ”ñ^½í½úå\‰ƒˆü›»½ï¸’âÛÅ]±M–‡­n¾¢ý„Ÿ¥‹¦œÃ‘[îÄ…½ò‰½½ö½î®¿N a“|çb½õ‹â¸…ÙêîÝ£¬Q¬n½ôð~šVæ¡üT½ðàtñÆñæŽ„„³ÚáàäÉ“ˆ²ø½û" },
	{ "jing",L"ƒ çRÛVˆiâ°Šn†™YÇGÂ€ìn„q­`îi¤•Ç¾®¾±Ï‚Š›H¾¹»~¾¢¾§Ço”ìŠù½›ìºãþÕeü ˆl¾¤¾©—}›Üû¶p‘ ŽÁìmÞŸÃ„¾­ö¦š„­Zåò›·Úê¾¦‚Š¾ª—Jã½¾¸„ÅÚåÙÓ¸„S¾·¾£•ßÝ¼¾´æº¶“ëÖ¸x½þ›¸‚åÉ¾°¾³¾µŠø·ÈªS‚ý¾ºŽyùX¸tìoó@¾¨ù~ØÙê€ëÂˆg«Eû—žDîK‚\ƒô¯džs“÷½U¾¥ G›G¾¶¾«¾¯„³œQ¾¬öLëæ­E¾²¾»Œc" },
	{ "jiong",L"½Ÿ•QŸ ƒ×Ñ•½Nìç‡åˆsïGåÄêÁã}ü†ŸKŸü¾½‚CÌSŸâ¾¼ÅQÌW›sîyŸ¡ŒlñoÞ›°›Óˆ·ƒTñ’ƒÕ EÅS" },
	{ "jiu",L"íƒH¾È“š™ã¾Î¾Å“[øFôbãÎ·c¾Í¼jðÕ÷Ý¾Æš”éN¼m¾ÌÈ\…Yð¯¾Äú˜ÍX¾Àà±¾¿¾ÃõíLèêôñ¾É‚w‘W`·T¾ÊöJÙÖ…Bšð ¬¾Çœ©¾¾®oèÑýnG…E–`û„óÅf¾Ë¾ÁŠe“A¾Â–ÍŽý›CÅi¼‘" },
	{ "ju",L"ñu‹JD ÊèLIé§™ÛßšêøˆR¾ÛÉXšÆ„ûÚªšjÚzùVÄ”¾áÞäéÙä“T¾ÑŸG›t¾Ýù‰’ºñÕƒh¾ÖÞê³µˆÏšÁÌ^õ¶ÅeÇÒ—ºÙÆÝ]à`ÈgÒz’¤ä|õáŸq„¡ï¸ÚkìÌ˜ˆ¿™Î¾åŸh“þº–‰±™h¾ãýeóM¶€‡ „—»Øe¾Ó×Þ±rñxèÛÍi•ZþŒøîÒ¾Ôè¢åáŒŠ¼‚eÔnÜÚ³^Ûg»Ûžì«œH¾æÁDïZ¾çôò‘§Œþ½Û¾Øø~ãIüŸ¾Ü—x¾×›®¾àÚ Šß ó†¯»‰àT›†¾Õ¾ÐƒâÂ`ùqö´¾Þ¾ß‚˜é°àu¾Ï¾Ù¾âÜì¾ÚþœÜM’±ÏJ¹_”Hâ é·ÎA’]¯Yé…ÜvŽe«~ÛBÝÏÂ‹„è¹ñÕ‡÷¶Ø‹õX„HÛRäzërÉ›àYŠèŠÛÄKŒÕ¸M‚IÈ{ñÀÉañ¾ÒÜ‡œ¦›ôúG’‡åðöÄÅ‰¾äõLÐýAÜÄlöÂ‘Ö" },
	{ "juan",L"×z®CÛ²ƒ¾íï…ùJÊ^–K’Ô‡üùNÁ\äm³ª™„æ±’ëhŽ†æŒîÃ¾ë—¨íj›û„»½çÄC…Û½vŸ]Çš€ˆ»ÑZðCìœ¾éÒN¾è¿xèðäg‘gÈ¦¾ìöÁŠ¤ä¸²CŠIÁIŽ™¾î„ÌŠFïÔˆ±ÈT¾êáúÃÄ–äŸïÃ" },
	{ "jue",L"Ó½^¾ð†­›‰…¨’¢ý™ùŠ¬œèö¾÷Ú‘ŒÖŽ@¯N™êÄ”Ïq™þþ‹²œÚÜ¾öâf“Þ½Å„äø_¾ôâ±ŒØÛÇ¾ó‘‰Ÿ]ñi u²Ÿš€¾ïÜjÒ™„]¾øÚbÚkšÜàÙÍXéÓø`¾ñÔEçå™@àåÄ_ÍDõêÐœŠxópç~õûŸØˆ« “ˆéQ‰øÓXŒH…É›Q{|Ø½~áÈìßïã¿›½ÀßIž½ÇÊ…’ÁØã÷Z‘ÝÜB«P¿”‚àÒæÞ³O¯‹Æ`”Çç«k‘•×HžŸŽD…ZÏpÞ§ØÊú€«iè‘¾ò¾õ­W" },
	{ "jun",L"ã—¬B´A¾þ¾ûùUðžÊ^¹êÐ‚ÍS—œ¿¤¿¥ðKûŠûŽ¿£ãzžFëh°˜¿¡¾ù÷•€ÞÜ›JÇq°—ÎD¾üñää]ÒŸ¹„…Íˆ­âx®¹‰ƒyùRŸóòE—T”|â¡õz”hŸaê}’löÁ÷å¿¢óÞ¾úÜŠ¾ýŠ®‘®å‹ùQý”ÈšŒ”" },
	{ "ka",L"Ð_ßÇ…í¿©¿¨Øûãl¿¦¿§†UëÌ" },
	{ "kai",L"†ËðÝažýžÍêGæbïa„’÷å|¿ªÌïÇ„ÑâýÔé_ÛîÝÜê]„P¿«âé‰NÐ_žG†þžÏæz™ü¿­ï´¿®•°ØÜ¿¬ç˜ŠKîø‡i" },
	{ "kan",L"‰Aî«¿°²™Ý®Ý‰{îƒýÛÉÝ¨¼÷¿²ÞRÐb€ƒÝêRíè¿³¸ƒê¬Ù©¿¯ãÛÝ|–Ý§šM‚°‰d¿±šK™‘´|ˆÉ¿´" },
	{ "kang",L"¿µ¿¸÷KØø¿¹˜±ÜãÊ¿»Ü{»~¿ºç_…H‡ã¿·³T¿¶é`â‚·^ß’‹¢o èýîÖ“•" },
	{ "kao",L"åê¿¼ÏË^÷Šèà›Ÿõw@êû¿¾Ÿ\”Ž·X˜‚õ‘¿½¿¿äDîíó}" },
	{ "ke",L"´Rï¾¸Sá³îWñ½äÛŽPîwÕn¿Å¯z“tˆÑ¿ÃîÝ¿Îæì ˜ÝV…Á„w³‚´žµL¿Â¿ÇºÇë´÷ÁçæÃmïýŒ¡ƒ¿¿Á…\ÚƒÄ´hò¤¿ÄåHƒ¾¿Àåíî§Á˜¿ÊŠÄºç¼ðâ¿È“U˜ÊQ˜}¿ÍËP”¨¿ÌžÜš£ä˜ÐŽƒÁã¡òò¿ËòSâŽ³`¾~¿Æš¤šÎéðáf„Ë¿ÉòÂšM„ÄÚÈd Éà¾" },
	{ "ken",L"‰¨ýlö¸Ø~¿Ï¿ÐñÌÃGØc‘©Ñy’õÃ\¿Ò³w«åo¿Ñ" },
	{ "keng",L"å”“¾„´ ¾³êlŠR³™äL’®ï¬ÕUçHŠs³nˆc“@¿Ó¿Ô" },
	{ "kong",L"³M³œì¿×Üwùy£›ïÙÅÁzåI¿ÖáÇ—¾óí¿ØˆÂ¿Õ" },
	{ "kou",L"„›“¿Û²]ADßµƒãŒtâ@¿Ú¿ÜóØíî·¸lÞ¢””„¼“¸ºpÜÒ²g¿ÙÊfúd" },
	{ "ku",L"Ü¥—üª@ýJ³L½f‡¿–ö¿Þç«¿áà·‚Vß ŽìŸ\¯‰õp–F‡ýÑF·”÷¼¿à³‚’H’¹¿â“‡Ñ¿Ý¿ßÚ¶sÃdØÚ¿ãÚœ" },
	{ "kua",L"¿ç¿èŠ¯•vÅ~m¿æ¿å½\ÕFîÙ¨ã’¿äóg†E" },
	{ "kuai",L"¿éƒ~•þ‰‘ñi÷ŽˆQôUä«Û¦¼[÷dX»á¿ë­g‰Kà””÷Ê‰Øáßà“ùëÚþc„SáöÒÄ’¿ìŽw¿êªœ" },
	{ "kuan",L"—p¸T¿îšLŒ’Œˆ÷ÅèwšE¸U¿íóy" },
	{ "kuang",L"…j¿ñ’[¹nÙL³m¿ôÜœ‘Çù\üY…NÜ’‘Èãk±q•ç ï²Ž½TÀkµVèkO¿õÝHÛÛD¿òäq„Áæþ¿ó½_à—û›¬¿ïÕNÚ÷ÝAÚ²›rßßÑb pÞÅ•pƒ—·ƒ¿öêÜÕE¿ð‰¿³qÚ¿" },
	{ "kui",L"—õ²Z‘èòñÌ€‹ÅŽu™œ•uŒºã´êÒÀ¡¿øÄCæKðŽhÛ“þ|ºˆ—óÖdŸåžà­¿÷óñ¿ùà°…Tþ}óY¿ûíŸòjÜiÂ˜ÙçšCîÌlØÑàkî¥Ëw¿ú‘|î`Ê‰˜æçqõÍÚóÄ„¿þÌwØ¸‹²z¿üÞñåÓÀ£»A´jÌðrè^ñùÝÞ áÂÀ¢¿ýÃvÂ‘šw‡]¸Qã¦í–¢„lš•êN…t" },
	{ "kun",L"À¤úAÁH½™òO‹Ÿã¶‘ˆÒˆÜÇ•‚Œ±åKéã§¶Ÿ³‰ÚÀ¥±—ï¿‹Gù{÷Õó‚—yé€öï‰×­@ÑTöH„õ«ÑXÀ§ãÍØ~çûÑ‚óˆµŒÀ¦ÎJðQØcîBÅC›ÙŠûdª^ÑhŸj" },
	{ "kuo",L"òÒœ‘²ìHŸŠ¶„îSÀ«ÊÊ’ˆ—IôUèéó–žN¹Q’•téŸípÀ©Èv”UíTÀ¨˜íAÀªÈu‡p" },
	{ "la",L" mö_ååÁÀ­ÅDÇ‰”jènÎ|´r“yíÇÀ®‡Äƒ•ËˆíBÀ¯ôFœ¼Ä—À²ñ®Î`Ïž–¬“Xê¹ÂäÀ°­†À¬—ïÞh÷véJ™Ê°]À±”YØÝðøß¡“Y" },
	{ "lai",L"ÈRª[”j‹@Ù‡áÁž|êãßF¬[À´ÕvÆ—…À³™Êù„†‹°]ÙlŽò‚g»[Àµ—®îsïªíòQžöDà[ån¹sñ®ˆœZáâ„Ðô¥²AäµÌDîmù`íù¹XäþÒs‚|üH" },
	{ "lan",L"™ÚÈŸŒG‰°è”ìµÀa”Ìž±íe¹ f»_Ë{‹ö¼hµf Lž°×Žè|À¸ÓEÀÃÀº”ˆÒ€ïCœ‹À¼­Šî½‰·°áYŽÓ™íá°ƒNˆh×EÓ[­sÒ[À¿ž™ïç›ÇžEÀ»äí ˆ™ìÜ_ A‘¾À¹ž‘…•ÌmÀÀÀÁÀ½ Šƒ‹Òw”G”rÀÂ•©‹û‡•‘ÐñÜÀ·À¶‡ÛÀÄê@ÌkÀ|ÒhŸ’À¾ ]é­ €»@" },
	{ "lang",L"˜¸”ïü†}„É†]‰iÕLˆ°ŸRÃž¹^Ý¹ÀÆÝõÀËÀÅÉvãÏæƒ–JÅ…òëÉ‡éšDàOï¶ÀÉäZ³„”ÉàHÀÇ‚Z¬˜Í™—O~ÀÈ–Tà¥‹™ÀÊò@Üq" },
	{ "lao",L"ÀÑõ²ÀÔÁÊ™Qñì›ÐÞL–U³ÀÒ·ðìç„ÀÐõu÷ÂgÇN³zÍŒ½jÀÎ»”°AÀÓßë‡ZªJÂSÜxèáÜ~‘Âç´‹áÀƒXã™ó€Âä˜÷†KÀÍîîºŒ„ºï©ÀÏ†ë‹ªÏo«™`“Æ„Ú†[ÀÌ‘Žî‘" },
	{ "le",L"˜SÁËö˜÷¦ß·«W˜· ¬Øìº{ÆIÀÕàÀßí‰’Aêbãî¸…³iàÏÀÖší" },
	{ "lei",L"×|Ú³°NèDèhïKÀ}‚ñ‰¾­zîž˜õªÀœ²Ìq˜ÃÀÙ™§Ì{‰CœIÀÕ÷mÙú½tË‰ÏœÀÜýF¶aéÛÞ[ìYÀnµˆ›æñçÀØÄB‰Í¯ƒ±ÀÚã™¦îL®šÀá›¤ûPÀßÀ×ÀÞÉ àÏµWÌrÀÛ´èˆçÐÕC…Ÿ™ïÀàÀhµX´ Ì…¿w”bÊuÀÝæÐåGî[" },
	{ "leng",L"Ýs±œÜ¨ËJÀä‚’´GÀâˆÙÛk¶ Àãã¶" },
	{ "li",L"Àåä‡÷kðÝ»šó»„îéöÛª™‚qÀ÷óÒæ²•åÞ]„˜”ÁÀôÑe³PÏ —˜™µÍjë`™ªû•b÷uŽ_Àùìc„°—ˆ”ƒžrÑŸ·ˆ‡ÑÀé ÓÌyÆnƒúÞ^ï®šs±LŒÞÏ[óööPØÖÉWíÂ™ÀÁ£†oçÊ»Œ‡³–^„^÷wµZåÎüŒVÀóøE•Ñøt”iµ`÷~Ï‹ôÏÀìÁ¦­€Ù³Àðž¦™æã‰ƒ«Î€¶]“à¬¶YòÃÁ¢”^¼H„{švæË¸{¬P×î¾°Oõ·¿rá­‰ÀõèÝÜV™ð÷óÉTäœÀü¼cæêö¨åGÞ¼ôfÇVèÀÀçûÀè÷¯–Û¾F—töâ—Ùµ’AØN…«ÑY iÀîžWárËžÜÂÁ¥ä‚Á¤õ”Àñß†ëx° ØÀøŒü²@ß¿ÀíÀúÁ¡ç\»hÝñØ‚ë_¶wƒ¢îºÝ°ãWŒCÀòÁ¨‡­ûZð¿ÀýÐGà¦Ç—°[Êkc˜»áB­–µ[ÀêòÛÁ§‡ÎÀï–ïÛÞÀ{…–”‰ìZÀöègóœÀæìåßŠÖ‚–säàækõŽ“…‰ÈúbÀûÀþùv›l sáû±XÍ‰WÚ\¹]­|Àë›ã•·¬—ŠÚðß›ÉõÈÏ~Æ´•²—ŸÅƒÎG—~«† À…Øª…“å¢‹KšÓóPÃš" },
	{ "lia",L"Á©‚z" },
	{ "lian",L"Á¬ŒD‹Õ‹t»^Ì_Éé¬æœó¹äòå€‡t‘Xì¡å¥Á°ÝüázÎ‹×`ïËO·SçöŸÈÌ`Á²Á®åbšaö–önÒœÁ¶…VÛšËW„ ™¹ÞÆ”¿ššŠY‹¼ÔœÇñÍ¿€›ËéçÁ·Á¯ºŸñÏÂ’Á­—†Ž R†ö˜öã”ßB´nÁ«Èjœ‹ÑžÁ´Ä˜À~­IôH‘Ùæ`“¢ˆäç Á³ÂIÁµÂ“Òc¬…Á±‘zÂŽƒI¾šÁª…UÖ‹»dà˜Á„iŸ’Âž‡" },
	{ "liang",L"Þc”Á¼Á½ÁÀåy†]†¤ìnÁ¸ÃžÝ¹Á¿ÁÂ†È’ëÝv†|Ü®ƒÉŸ´÷Ë˜ÅIôu¼Z›ö‚ZÝˆÁ»é£¾HÝgÁÁ¾nÙûÍ™ö¦‚zÁ©Á¹õÔ‚ŠÁ¾ò@Õœ´ÑoÁºÎW" },
	{ "liao",L"ósÁÊà€•ÅàÚéRÙ’ç‚Þ¤‘’ºƒ¿ÁÎÁÉ rÁÍÁÃÁËŒ×ÛŽË€ÁÌ­V\Ä‚¯ŸØIŒ®ÁÏ‘l˜Í²tÁÅÁÄÞÍùß| vÁÆâ²çÔïmïfúÏiÜGîÉ‹»ðÓúžÒå¼á‘¸NÏYÁÇŒ³”¶ÏoÁÈ¸XÛÄkƒJx" },
	{ "lie",L" Ú mÁÑ’ž±ŸŠGø•ÞæÛø¾FŽ{ªdä£«Cƒ•„ÃßÖ—˜Þ˜š¸ iõñŸI÷àh†`ÁÒÙýŽ_õhŠ² MÁÐ™§Ã‡ÁÓþš’£Ÿ­ÁÔ›¼è÷vïVˆ´ôóÂ~”YÆ”ôQÍ}" },
	{ "lin",L"Åõï[ò•ÁÞëO•Éƒäà®V¾D¿šâÞçlì¢¹ƒÁÖ™_¯rÞ`ÅR‘¬ÁÙžŠÜk÷[ïCÁà÷ëõCÌAŸi°»‘ÞOÜ\ê¥”Ý…éÝÝþtÁÛ…›åà‰É°SÁÚ‡ÁØÜCª­UÁÕî¬á×ÁÜôÔ„C°R®ž˜ð«ãÁÁ×ßøÂLÁß´@“ÔÈHÙUéŠŸûƒj•—·Azû‹ÁÝ" },
	{ "ling",L"ŽX‹ø™ÐâÊ™ÅzûwÁí’èžƒû™Áéˆ{¬O‰çŠêê²’°s•`öNàòÐeÐ‡Ôfñ|êtõCÁë™ôèÚû_žâÁç ‹öìÚš„cîIëžë‘¶{ÁìÁèý’ÁêÀâëãöýgáÁã¸ Éˆ ÷éqÁæ©–Š–ÁäøoÁáôáì_œRë™¶ èù¾cä™ñö U‚’ÎÁå¶ÑkÁâ–EÝsßÊÝCÌhHµ’¸nýhÁî³gÛ¹½@òÈÜßç±ÊCì`" },
	{ "liu",L"´e–Î¬–›fÁò°@ìÖ”éŸÞò˜ÁóåÞ®‘‡®”åöÌ®œëw¾^´zðÒÁñ®qÁöðs‘Ëúwç¸ïvÛ‰çsñ‡«€æòäìCö†—BÁøÂµôjï³ÁõçBæyÉsÍÉ]òtÁ[AÁùûméH˜ñÁ’™PÑ^„¢­]ÁS—Pïiïfê‘Áï¸Ë˜‹ˆì¼žgûˆÏYïÖ‰gñœÁðúV‹ô¬ŠÁôä¯Á÷Â½ïd" },
	{ "long",L"ÁüýVûTýˆ¸_µb”n‡µÁýÁú±€ØF‰Æº\˜™ÐH—YÅª zóGÂ¢»\‰ÅñªÜÂ£ÐFÃ@Œâ³ŠãñÂ¡•îŒ™­‡ƒ¥ççÁûë]ÚL–VèÐ†UÜ×¸oíÃ™ÉŽaÜ[ðèxÛâýŽµaž{ŽbÂ¤ì_ØLìNçXÐiœ¬ÒtëÊ•o²”ÁþÌdÊ”œö" },
	{ "lou",L"‡D¯›ÏNà¶ºtíVïÎÖŒ²kÂ¨þoUÂ©âÊVÂªótÂ¶ÝäÄ|÷Ã¯œñï“§®RÂe‰vŒÍ¸Mœ¾ñÀ˜ÇÜ}çUò÷ŸÓÂ§ƒEÅ”ðü‘fáÐŠäIßsÂ¥Â¦ÙÍ" },
	{ "lu",L"œGÂ¹µ“®f˜ÌåjžoÂ¸­ˆè´”]åÖ·c„Î‡£éÖß£‰ÀˆP—¶Ì”Â±Â¼ÊIFûRÊ€ t²FüuÇŠ™¾ä›Å›†ë»œ’ÇÌJ˜Ä³”÷wçeëªÂº±RðØ„—º—çG ê‘ÓtÂ¶‰oÉÀrëÍÝ`´{¾GÏFÂ­èuú˜»U²’ŠáÄyÖ}žZÚ€”d–›èÓ“ïô”‡´™©ïBààÆ@»VÂÌÂ´¾vóüùc«GÞ¤VÀžéñïåLùnÙTòƒˆvÂèz³tÂ»Þ_¹‚æ”ÛäºŽäË«SÂµöIãòåhÛÂ«Â¯Â¾ÁùÂ³œOÛjéûÁ’äõ‰nÌŸÑÂ²šÚº˜”mÆAòJÂ¬Åy÷|Â®ô—à~Â·ðµÂ°ê¤ÄrôÞA¶˜“¦­oÄwózÐBöÔÅF¬fôµ’ çœûuáX]„ÛÓ€±Jøšâ„ƒJÂ½" },
	{ "luan",L"Œ\ž¤ÂÒŒDÁcÅMˆKÂÎ°f™è°gËHð½ÂÑöÇÂÍ”æ®ÓTž´•ðÙõùFèïÌ‰ÂÐÃ‡ŠaÅLèŽû[ˆJyá›vŽnÂÏ" },
	{ "lun",L"Âbö´KöMêÂÖ‡÷—‹‚Ý†Ä@ÎFÂÙ¾]¶—´ˆÂÛ¥ÂÕ‘ÛiÂØˆÀ’àÂ×Õ“ÂÚ‹EàðÇ’†íœÓœSä—" },
	{ "luo",L"ÞûÌ}«MžTòŸ„sÎÂßÓZé¡ÂÞÜs½j”{ÓTß‰ÉzR’Ó•ïäð‡ÓíÑ q³ŠÎÏ° ÎÙÀ°eãø¹JçóöÃÂÜÂáÂãñ§Üýµ[úŸèŒƒ±”‰Âæ†ªâ¤ëáÆŒùBÅIƒ¬ÂÝæ ò…æƒÀÓÍxÂåñ˜ÀzÄT™åÙùÂç™ïïÝÂä³`ÜVÂâõÈõið”»j› îbÁ_ãtÞÛÂà" },
	{ "lv",L"¹˜•ì·tÂÈñe¿†ŒÒˆ‡ÂÌ¾vÊVÂÀÂ¿Ÿf†`ÂÁé‚úyƒ–…ÎéµÂÆžVëöÜ}¯È„…iàLËƒÂÂ‘]‘fÂÅïùÙÍóHäXñÚµ~Ò@ÂË l™¬·„ãÌ¿|Ä|—oÄoèr¾GÂÊ™°Œœ½…ÂÉÂÇƒEšÑ‚H„íŒŠÞÛÂÃÂÄ" },
	{ "lve",L"ï²äx®ˆÂÔ·DÔ›ŒœÂÓ„…ˆGäs”^" },
	{ "ma",L"Î›ÂîÂðÁRßj Ðó¡‚ØñˆÂìúiÄ¦æ‹ö‡²K¬”áïœÔ‹°‹ŒñR¯rè¿ôK˜qŒI¯qÂê†áµlÊhÂèÂéÂïÄ¨µTÂëéUÂí¶MÏWªw´aò‡" },
	{ "mai",L"Âò‡OÐ]ÂñûœÝ¤ÙI„êÊ{Âöú”ß~Ï‡XßéÂõØ‚ÂóÙuö²Ã}ì@ÂôìAÛ½Ëh" },
	{ "man",L"ÂùÂñò©ÙÏ\Âþ±”ì×÷©Ü¬à„î”ŽÂýÂüƒK˜ÑÂûôMŒÌÏŠ›œº˜´ÊAçNÂú‘`òý²–á£MÂøÏTæž“¶çÏïÜÃ¡Ö™ÐUôN÷´ðzö ‹ ²mªƒ¿zÂ÷" },
	{ "mang",L"ýýˆ®mèšÃ§Ã¤‰Üšû±ZíËäÝ–M……¸ˆ…¹âIÃ¥Ã£³‰Ží–xä€}Ã¢ªKñ Ã¦ ¯Í{Çƒòþ†WÏ‘ÆŸ¯gŠÁ–nŒ´ ½›ÀÌMÚø" },
	{ "mao",L"Ædó¹FØ¶Èr‘ùÁEêÄêóÃ¬Øˆ‰îè£à|¾ˆë£Ã²ÎcÃ³ ÓÃ­•§°páFÉ‹ƒÓå^÷ÖšÊÃ¯ãTÙQ±g‹uã÷í®¶mÃ¨á¹Ã°—ûàŽ–‰š»Ã©Ã®ƒØˆéœ~ÙóÒ‘î¦ƒÐÃ±ùšÜšó±ìWÃ«Ãªòúãwì¸šÓÜâ" },
	{ "me",L"ŽÛ÷áÃ´üN‡ª" },
	{ "mei",L"…Ðœ„ÃÀ¬CàdÁoÖišîÃÁÄPÃ¸›i²‚ïÑáÒÃ¶¹Ÿâ­ÃÂ’¯²SÆ€Ã¿íi¬sš°±ÃŠ˜ŽŸ¢œŽÛÃÃÝ®–Ï›]Ã»ÃzÌj‹‰æ[Ã¼úBÚ›±gäØ÷ÈômñÇÃÓ’{¶C±tÃ½Ãµé¹¯cä¼üqÃ¾±ŒÃ·˜M‹ZæVÉBüeˆõ™­ÃºˆbœÕÃÕƒñµ|‹Ê‰r”uðÌÃ¹ BÃÄÎnäY" },
	{ "men",L"•¹ŽÞÑå{ìËîÍ’Ð‘¿«j­JéTÃÆBÌŠÃÅM·`š‰ Fí¯žÇ–éYÃÇ‚ƒ" },
	{ "meng",L"ÊpÎ{ÃÍÃÏûLãÂìXŽíüw’úòìõ’‰ô˜ý²“ÃÊëœ‹“Þ«ÛÂ‘ºÃÉëü«Bð•äÇm®mÝùë‰òµÃÈô¿ƒáöQåiíæƒûsü€÷‘¸Ã¥ìDà–ŽÌœÉà‘²‰ÃË”BîŸ®HìWÃÎŒ´ô»ó·ÃÌ‰õšÙ" },
	{ "mi",L"åô÷ãûJœPÊZôÍÃÑÃÛº€ü†™Ž¶Éo †åµôéÚ×Ò“e‰Q·`ÞÂŒBû†ƒßŒªáˆìòÃÒÃÜ¶[ÉqÃÖóÃÚŒ©›^Ök˜ÆµzÒš‘ÛDœ}”}±‰»…²yÁ]Öi÷çÃÔÒ’ÁdÃØŒsëßà×±~ –²[éSÔ™­Œ“º˜až…ÈŽÃÝ†OÊUèfŽÈãè›mÃÐ›aßäÃÓƒçá‚«J„¯ãñ™—â¨ž§ŸÇû”ðÃÙÃÕØÂÃ×áƒ›¦Ëz”C" },
	{ "mian",L"…›DÃå²õ| ¤‹išóÆÎeÈxüIÅXÃàŠå†»ü@¾‡ö¼äÅÃæ‚aìtìr™†‹î¼E–uÒÃäüwÃãÏŸ‚Áû çÅ¾‚ÃÞÆP²Šëïíí¾düMãæ™¡Ãâ½ƒÃßäÏ¾’²Œ„ÒÃá" },
	{ "miao",L"÷]ÃéðÅÃê«Qç¿¿ŠåãÃçÔN¾ˆ³®ÃíŽøRÃìíð¾˜£¸kÃë‹bÃèçÑß÷èÂù‘ânÃîºFíµ‹·" },
	{ "mie",L"øpßã…¸ŽÏžf†_óúêEØ¿Ò”÷xóº™­žûÃðŒP‘Ì±uœç“}Ð`ËIèfÃï" },
	{ "min",L"´CÃò’Ï‰‘O•Gä æFB¾‡áºö¼•¡ÃôÁFœbëíªÃñüw‘‘¸œ˜éhÏŸé}„ÇÜåÃö¬zøsƒí„bãÉº‡Ãóçë¬Yˆ„öšçäœ¡³R÷ª”•œ“çÅÙ‚…Ý¬\¾rƒoŠ“”°•FÃõ¯xâŒãý±a" },
	{ "ming",L"ÔšÚ¤ÜøªuÃúäéÃùÃ÷Éqàp–Lƒüã‘Š±±…ÓK˜iøQÃüÃøõ¤êÔ›³‹“±b‘DâÃûî¨" },
	{ "miu",L"ÃýçÑÖ‡¿Š" },
	{ "mo",L"÷áðxÑJõöï÷âÉéâÖƒª…ºÑ–£šzæÖÄ°±‹Ã}Ž’ÜÔ•b”VŠ‹õøŸo»Š jñòÄ£\ôŽÖ„Ä§×O•½²aæÆˆ\±uÄ¢ÄªïÒÃ»›]‰s˜í¿}½]°t“á½QÄ¯ºÙÄ­ÎÞÍˆÄ¤Ž”ÅžfÄ¡ÍòÃ°²hÄ¦„¯ÂöüaÏ_Œ­ñ¢„¹æŸ‡±Ãþ¼UËüNórÝëÄ©ƒØÄ¬Ä«Ä¥Ä®ªC†ùÄ¨³]‘ÛØ{š{ÀgÚÓ‹º‹ßµcØ€ã€ì…”}‡¶‘½Çeò‡ð‘±‰" },
	{ "mou",L"²y…ÞÏwÄ³íJ¿ŠcÄ²Û_Ö\ˆéÄ±„À”–øœ›£ßèíøçÑüEòÖöÊ–üÙ°ãw" },
	{ "mu",L"ÀÑÄÁŠÃ‡`—úÍ]…Þ ¸ÜÙ®rÈrÇ€ÅÄÂíJžÑØïÄ·ãåîâ®Ä²š»ë¤®yÄºÄ¼šÒŽ¿Ä¸ÄµÛéÄ¿Ä´ãfÄ»Ä¾Ä¶Ä¹˜ÒÄ£®€ ñ–]‰®ŽÄ½ÄÀ³cãaëŽ‘Hë‚ÆŸß¼\¿}Û[Ãk" },
	{ "na",L"Š{ØyñÄì„ÄÆÄÇÔië~ÉiƒÈÜ˜âcÄÃ’…ÈÄÄæ“ÄÅ’‚¼{ÔFÞàÕyÄÉ†òô›ÄÏp¶gëÇïÕ¸™ÄÈÉSØvÐœ" },
	{ "nai",L"ÄÎÄÍÙ¦Ø¾á‚™ÄÄ‹èÄËèÍÄÊÎ—Š…ÝÁÜµÄÌ¯Gœ‡ÂYÑ”åriÞ•" },
	{ "nan",L"®~ÈlŠ{–•¨ÄÐàïòï–¹Ç~Ÿ²‹R‚Oôö–ŠÉÖQà«ëîÄÑÄÏœ¯‹©ëy“D‘Ú”‚’oéª" },
	{ "nang",L"“rÐL™ò›ïâÎƒ²“îôTàìÄÒß­êÙýQ~ð–ež²‡°" },
	{ "nao",L"«LÎjFÄÕŽH‘«DòÍÛñ´Zô[“ÏçtîóÀÔiâ®Žjˆß‹š´LßÎÏuØ«ÄÖíÐ‹CÄ×émÄÔÄžÄÓÄXp…D‰ë×Dè§" },
	{ "ne",L"ÄÄ±„ÄÅÔGÄØ’fÚ«ÄÇ" },
	{ "nei",L"ŠÌšß›ÔÄÙð]Ã•ÄÄðHÄFÄÇÄÚƒÈåMõƒõ" },
	{ "nen",L"‹\ÄQ‹¯ÄÛí¥" },
	{ "neng",L"ÄÜ¸o" },
	{ "ni",L"Ää’vîêÍ‰Þ‹•ÛC¶vÄÝÝrœNñDËoÛèšî’í¯[—´ÄØÄß’f™ŒÉ‹¤ž—»uýuˆÐ©ÄâëWÄàŒÛŒN‹òÕyŒTêÇÙ£ƒ¹ì»ÄÞãbÍeÄáâ¥Äåöò›mÅMâõƒºíþŠöØƒÓrà\¶[ƒ“èXÄãÄæÆsÄçâ‰ÄÄQÂžÃf–«ð•¿öFˆÓî¿QŠ…Ã”M±z‘¹ ùƒŒÎUûŒ" },
	{ "nian",L"öÓŠ¨õRð¤Ÿˆ…`¶jÜTÚ™‰|†P“ÓÝšöTÕ·Ý‚Ûþ†ˆöóÄíØ¥›ÝœVÄèÛœÄëéýÅˆÄî”fÄé¶|Õ³ÄìÄêºv" },
	{ "niang",L"‹úá|ÄïÄðá„" },
	{ "niao",L"Äñ˜ÒÑ™ëåæÕÊ\ÜàŒ³ÄòÆ›øBôÁÄç‹Ø‹–ÑU" },
	{ "nie",L"ÛWÄõÇŒ†ÇÅYô«¯[—´äŽ“µ¼fÄù™ÇäO–¨ým“Iœ¸Äøæ‡ºQÐAéÞãbˆ[à¿¼b˜®Ûh‡ËŒZŽ‹ïDêŸèX‡Üò¨ÜbØ¿Äó»HÆ}è‡ÙõæÞÁ”¤¶Žqè‡ÄöŽLÂ™×‘ÔÛfåR‡§ÄôÚíÄ÷T" },
	{ "nin",L"’ŒÃ€‡áí¥Äú" },
	{ "ning",L"ÄþÄüÅ¢è_Å¡ôËf‚ž”Q‡“ßÌŒ‚²…ÃÄûñ÷ƒ‘ØúªŸ‚Aå¸ÄýôXÆrŒŽûH™ŽŒ|‹ÞŒ„Œ‰ôVÂœ™F" },
	{ "niu",L"Å£›S «FžÈ–ƒ¼~ìÇyâoâîÅ¥áðÅ¤Å¦æ¤ÞÖ" },
	{ "nong",L"·výPÅ©’° \Ê¶ŒÒa™×ßæÅªÅ¨Ç_‡’˜°Jâƒz™`Ä“ÞsÞrÅ§áxÙ¯¶Z" },
	{ "nu",L"ñw³eÅ¬Å«Ôi‚Õ¹“xæÀÅ­æÛåó¹@Âæå" },
	{ "nuan",L"œ¨Ÿð`Å¯ŠfŸœœq" },
	{ "nuo",L"ÙÐÖZ“—jÞßS¼K’ý™DÅµÂXÅ²ˆëšÃ“xÛ¼XàG·zÑD˜`ï»Å³Þù·LßöÄÈ‘Âƒ®åŸŽ]Å´" },
	{ "nv",L"Å®ô¬›\í¤âSÂxîÏ»s–HÐZ" },
	{ "nve",L"¯‘‹FÅ°³–Å±" },
	{ "o",L"‡—àÞÅ¶" },
	{ "ou",L"Å¹É’®T™¯‰p”·…ËÅºšªøkáq…^ŸàÄUÅ»šWaÄpÊq‘YÖŽý{ËšÅ¸æ–útÇø‡IÅ·…¾Å¼Ïñîâæ¹pÅ½Ú©ê±" },
	{ "pa",L"èË°q°ÇÅ¿Ð’ÅÃâZÆtóáÅÀŠr’öÅÂšñþxîÙÅ¾ŽÚ•ÅÁÝâÅuÅÉ°Ò" },
	{ "pai",L"ÅÇÝ‡ÅÆæW¹uÝå—“ÅÈº’Ù½ÆÈÅÅÅÄªT ÛßÉº”ÅÉœkßß" },
	{ "pan",L"ãÝ±eŽ´Ég¿T–®È_˜„ZîGÅÍ®‰ÞÕžÎ·¬ñá´B°âÅÌÎŒÛA žãú“„œ°íQÅÐ´‘±~ÅÖÄ‡ÅÑÔj‹ŠˆmäƒñÈÅËè‹œãÛ˜›cÑ—žbõçó´°ãÅÏÅÎÅÊƒë›æo±P" },
	{ "pang",L"…€ÓI›PáÝ†çë„Ãp‰â÷›óo‹˜°õåÌäèºUý‹ó¦Ú“Å}ý‰ÅÖÅÔÅÕÅÓ°÷ÏìQö„žÐÄtÝò°ò›`ÃT ¥ÅÒæ^" },
	{ "pao",L"ÅÙ­”ûƒüB‡¥·•µPžäÅÝÝNëãâÒË‘è˜°’Å×ìŽ·…ÅÜïRðåˆƒÖcÈaÅÚ’û³h ÅØáóÐˆãE ÜŠEÅÛµ^ÞË" },
	{ "pei",L"àÎ«˜ÅæêkC¬ašÅÙr’yÅá›ÖÅßö¬”ä éÃSäžÆž¬‚_ÅãäÄÞ\Ð[Åäñ]ì·ÑpàúŠ³•^ÅÞïÂõ¬ÅåïÅàŠçÅâ" },
	{ "pen",L"Åè­›äÔÈ†Ðv‡ŠÅçå†Ïš\" },
	{ "peng",L"œAÅî³yÅíÆMymÜ¡ñs’üÅò“žó²‹íŠÝ~íŽàØ´yó—¸†åAâñéoßJ„™—ZÃg—Õ’²ÀeÅóÝJÅëÅìóŸÅðÅõ‰kòuŸÔºU±ÅÅé¯nÛsùiÅêÅñôJp‘u˜¨Åï°vœKÇlÏe„ú·@›€“s‰X˜Õ—ÄÅôÌX‚‡Åöèm" },
	{ "pi",L"ñâÄmš·’y‘š²DŠvØò¶uºfÅûãÚðÝ·K÷‰Ü±ãYWµFúÁ`‡‹œŸæÇô“ó‹Û¯î¢þžÆ£ÚüÃ˜ ø ò»µñyêoØuò·ÁTÆ¢õQêVê¶š¯wèÁÅý·ñÆkÕ|äÄ¹vBÎ“”è°Õ±ÙÂ\åCÃY¼„‡º¸“»zÆ¥¶y‡š³žÌõBî¼âWÛýÆ¦âtÜÅÍnÅøõùÏKÄMøaÛÜñÔÆ§âÍoªWÐK‡›Æ¤´iÁ‘–C—Àòçã›ÅüßÁïâbØwàèîë–ŠæqÅúÅ÷ÅùÅþäšÆ¡âÏØ§ç¢‰ªÆ¨µGñ±“FÚéëRÉÆ©¶ß¨â”¯@" },
	{ "pian",L"Õ—ÒÆ­ÛMÕ›®Kñ‰ÙGòNÆª¾œÙX‡æú@ò]çÂõä­p‹xÞq±ãò_Æ¬ôæåw˜FªpÚÒêúÆ«—èæé±âÄAójëÝ" },
	{ "piao",L"ÊEºgêQÆ®áoóªÆÓïhéèØâôwî’„Ü¿~ƒGàÑ”ô®óTÒïg°ŽæÎ‘GÆ±Ëi÷ÔÝ³æôÂHî©Ö€ ÜÆ¯òŠÆ°óQçÎ" },
	{ "pie",L"ë­Æ²“ÅÒ”•ÈçvØ¯‹±Æ³‡ÜÖ" },
	{ "pin",L"·|óDæÉŠÐÆ´Æ¸Æ·ñPÞÕµI²‹‡¹¬VÏ™šý«nÌOò­æ°³WØšîlËdêòïAÆ¶‹åÆµé¯" },
	{ "ping",L"³yõG‘{ÍƒÆ¿Íg¹’Æ½Â†«rÀùd™qŽ—É‘ÆE®JŽ±ŒÎãuÙ·ñTîZÃg›¯…çÆ»·ëèÒÆÁ›ÚÝƒàZæ³ÔuŸvÆ¾ÆºÌOÝZœKÆ¹³f¸z„R‘k›€®jŽ£öÒÆ¼ºqÇLÆÀ" },
	{ "po",L"óÍÆÂñpžTÆÓÆÃ°~ÆÉ³kÉbŒžáNŽˆþ‡ÆÅ“„ÊXçêÛ¶ŠŠUãOÆÈÆÆ•^ãøØÏîÞñFœ_ÆÇœ”g™†\FîÇ·±îHŒûçkê·ªtÖc”’›¨ómŸBÚéáwÆÄá•²´ð«‡M" },
	{ "pou",L"ˆ¡’gî_ÆÊ…Ä—”†V’½à^¹r ÁÒJŠËÙöÞåŠç…ð" },
	{ "pu",L"ÆÒ™kˆOàÛë«“ò•®ÆÓÆÏ±©ÆÑïäÆØÇ[Àb’ÃÖEÆÔŸMŽ}ÆÌÙŸÆÕÒiä„ƒÍ—¸¬×VªŽ‰áTÉhÆÙÙéÆÎÒLè±ïèÆ×˜ã–¿¯jÆÖƒWå§‡þ²rÆÍžÊ±¤·o“äçhõëÅmäßç’õ‹ÆË°þÆÐê†ª" },
	{ "qi",L"Æò‘sËjŒóì÷Ûp…æç÷–ÖÆî”ÆŠí†ƒÆïçùÎB»üœDí¬ÆýÜÎÆúÆûùušâ°žÚ|››ë’Ù¹™–ÃX÷’òà‘hôGôtö’Ÿdƒ[´wŠÝÇKòUžÅ™ûÜeÕƒÆæõšÍ[èŸÆãÆôÍVÆø…Ñá¨ü±[†¢œæÛe’åõlÆ÷ì¥ÆüæÆÞÆÝòÓôë«Où³žÜùÍTÆßÆÚÆèýRÆêÌIçˆÜjÆóõèÆZÆçãàù†èç¾N”Œ˜Ú–¾e”çâHù}†šÜ»ßŒ‰ó¸gàVÆáí ÷¢¢­DÆñ®P¼©Ü™ªXÆðñý‘iÝÝÖHÆõÑv•´÷èÆå†uˆÎ¾_³H‚ˆÏ„ÑEÀdè½û˜àÒ—tÆéôn…ý’Máª´ƒÆíÎ‰Ã´JÆëÜ•»žÅ Ï“”ÅŽÜ—¤“ ÆÜ•’ØÁíÓòTÆÛÅpô—R—‰¶QÆäÆì¶ò€Æâç²åWôy…ÛßçK„~æëä¿—Ž’ÝµJÃIËsÄš«^‘¼êÈÓsÔ¯O—«Ý½™‡ØMéÊÞ€ÈWôìÞ­Ñw”ª¼–†™êMôoÓ™‡rÝÂ‚úàœÛaµoœŒÆù´\Æö´„¾LÆà¾ƒ–OÀ™" },
	{ "qia",L"³L“üˆXƒî÷Äš“U´l’‰Ç¢ÝÖ˜HáM³sñÊgÙƒr¿¨Ž˜ÚžÆþÇ¡" },
	{ "qian",L"ã@ÊgËübÉ`„X‚¡Ç§ã»¿yâT ¿ò¯ÈœÇMòcqžÇ´Ç¶ÎS¾P’®‹ìÇ­åXˆU‚ßƒLÒ¹ˆìy‰µã¥ßwç×ëÉÍO”pQŸÈÇ¬½šþ™NÇ@•üóéØ@øZŠúƒŽ…•á©Ü‰qå¹ºR“îväEÈ“‰‰–e‹üŸšškŒòºž“bâjÀ`˜ Ç²Ç«Åˆ”oÇ¸ôSÁ{ÊnÇ¯èýÜç™Œ†k”q°|ía’R÷œÕ¸dÇ£çc‚]˜pûeãUôRòqŒÜ·úY×l RÇ·ø¤ÇµëeêùÝ€ÇªÇ±ãQÄd“Ã‘aÜÍºGœ\Ç°öÃÇ¨ÏËí©å½›F™÷—Úä“¾Ç®ŠdÇ¥’ƒ™¥Ç¦žKÝ¡ö‘è~åºùkÛÉÍZÙ»Ç³’Š˜ÙÝåDŒRÇ¤Ç©ècèBîÔ‹`æPˆTÅO’ç»`Þç’Lâ`Ötò`" },
	{ "qiang",L"ª}‘êÇ¿û]‹ÔÁmÛ„èÁzZ™—¾ ›çIãÞñßÀH¬šÇºŒ¢†…‰¦ìÁ¿‹Ç¾“°Ê@«oçjÅšÇ¼¬j\ÁuäÛ–ŠœÙôÇ†óŸÍª]ïºïêÁ†ê¨ËNÌbÇ¹“Œ˜Œ „ßòÞ“¬½«ÇÀ†ÜúIõÄæjïÏ‰‚éÉÇ½”ÖÇ»æÍ™{Öm" },
	{ "qiao",L"‰”‰U“³¸[ÇÄéÔ´íÍÇÇæ@‡a½¶å ˜òÇJ¹›á Ú½èAÇÁã¸î˜õÎÇÍàb„äîNÕVš¨ófã¾Úˆ‚¸ó~ÇÃÜFÊwà…íXŽÉÇÅçØÕÐ“ê±š¤ÇÆ´`‹´´“Ÿ÷™]ÇŸÇÈímçy†ÌÇÊó|àƒÏf¿ÇØä³~‰ŒÇÂþŽ»ÚÛ¿”ÂNá½ÇÌ×S‡„ÜNÜE¯ È¸Ÿ}Û^÷³ŸòÇÉÚ‰ƒsÇÏš£ÀRÇËÏË–Üñê~ÇÎíIŽ‰§´™ƒS´x" },
	{ "qie",L"Ë~†éÇÒÔˆ¾fóæôòã»ö@Ù¤›­ã«çƒôŠÇÑ‰êüÇÐ–A·lŽ¨Ü‚Œ…‚öl…L›ùÛoæª¸›ºDïÆÍ„Û§¸`ÂÇÓÆö¯CÇÔå›°m" },
	{ "qin",L"Œ€û²›ì€ùj•T’aŒ˜†wôÀÈBÂlÇÞôgä“lÏˆÞìâ†ÏOÜËÇÜ”ÜÂ™NÇÖõˆ¨ÓHÇÝÝÀ¬l…ÂÇÙäuî›‹]‹Žï·â\ñŸÇßßÄôÌC¯‰ƒ’RÚ_ÚcâdÇ›šJâÛƒ¡Ç×¸øòûÇÚ‘¦«ˆaÇÛàºÇ™žpéÕÕW‘¥“å‘[øVŒ‹…ñûÇÕäÚ’ÍÇØîzÍZà½óVëdâs—v™Âˆ²àß" },
	{ "qing",L"ƒ Fõ›ÇêíàÇáš Ç×N‘cœ‚¾PƒAÜÜœ[ÕˆÇãÇëÇçìmÇì•¦†¦š„Žö‰ð³ àW˜½ÝXÇæÇmþ‚—³„ÍòßÓHŒx÷ôÝpÇà„…ƒõþƒàõÈéÑ®_è[óäÇåôìÇè™”žDìišäÇäí•Çé“÷öëö¥ÇâóÀ’áù‚ˆ½" },
	{ "qiong",L"™KÇîŸÅŸz…o¸\ÄŸ¦Ÿw­Ž¹HöÆŒ^Çí­‚Ë}–÷óÌÅ|²`¸FÜäÚ^ñ·þ{Úö¬IÍ‹õ¼‹Ö‘wòË±žË•ƒ’­W" },
	{ "qiu",L"Çô½‡òø­GÙ´öq¹êäÐÈcÐ@ù”ÙgûjÞíFŸª÷GÍÇiÃF’ºþ•Çð·hnÇòÇóåÏÉ’«U†pºE¦¾ÓaœrÇï¼zé±öú–_›½åÙ’@‹pšÂ…´ÇñêäÓpáì™Ïöp…œÇõ÷Aˆw³ð°“šüÛÏÏbäM÷üHÌU˜þÎ~±HòÇÚ‚Ó‰šðý•—W¶k“zÓˆá ³Ž€õ‰ÇöõFábíGò°ôÃÍAš‚ôÜý”á–œª" },
	{ "qu",L"ž›êrÒÚm»cÂ^é˜ÐR‘tIüšèLŒþ…JíáùŠáé…é‰…^¼ çÈ¡’|»–ýxð¶Ú ßÔx‘óàTÃlõ@ÐdÁ”Ç÷ÇþÎƒìîÇüÚ…ñnÇ­SÜdÃaÚ¹LÇý½Mó½ÏJüLÈ¤ó”ÕoÝ@öÔsãÖÐçEÓYˆoÛBÛ¾Þ¡È¢”×Þ¾…í÷OÅJ”·ôðÇùÐ ™áñ³ëÔ„`òŒÈ¥›µŸaøzÚ°ÂJÓNÇúÇûñlœTüD¸yèŠ…Z½P÷ñë¬Ü|È£ÓUÏgÇ†ÇøüzÍmöÄá«¸lè³òÐüCõLÜÄêïlûYÈ" },
	{ "quan",L"®lÈ§ç¹†­œ²ÛIçzÈ¯ŠºïEÝböeéúƒZñî°È¬›§™à÷™‡ü¿X„ñÜõ³ãª¬†çÓjýjà —¨òg »íjÈ¨Žk„áÈ«êBÄCÈ©È°˜ØÈ®Ì†ˆ»È­îý½h ôžµ“‘È¦óÜ÷Ü ºÈªŠIŠ÷¾J²•ÐSw³oÈ››LŒA ÅÔÚ¹ãŒžïÛmòé" },
	{ "que",L"âÉUª“‰”µC‰Uš¨ˆ«È´È±Åb‚àÈ¸”¦‘UÈ¶é í¨µ]Ú|Èµ…s³‚ÍXÈ³êI´_È²ùo P°”“nÈ·ëa‚´`ã×´FãÚ" },
	{ "qun",L"n÷åÁt‡ï¹„ûŠšVÈºÝlûŽŒlÑdåÒÈ¹‰æŽ " },
	{ "ran",L"‡YŠ˜‹vó†ÅjÈ¾È½™LÐ€÷×Ð…È»…ß¿‘òÅÈ¼–¹ÛœƒÑÃV«zÐ™ÜÛŸßÍc" },
	{ "rang",L"‹úè‚Ïâƒ¨ÐLìüôX«KžÈÂÈÁ×j‰´ð¦ÌZ }„ðÜ`×Œ‘Ó·yÈÀÈÃÀvÈ¿" },
	{ "rao",L"ðˆëNÒYÏuÈÆæ¬ßvÈÄÀ@ áÜé˜ïèã”_Ê‹ÆÈÅ" },
	{ "re",L"ßöÈÈÈôÈÇŸá" },
	{ "ren",L"¼Œ¶‰¶eøžÓ•ÈÍÈÑôŠž ®ÕJÃMñÅïƒÀÑGÆ\ÈÐáÈÊ–ZÇŒïš—eÈËâJâmÈÒ½V’PìzÄHÈÏã…ïþ¡—ª›ÝÈÎÈÉÇYâ¿µs„Uì~ØðéíÜó–ß–kÜ×šŒã–á¼xäÝØígÈÌ" },
	{ "reng",L"ÈÓÈÔµiÆeÞwê—" },
	{ "ri",L"óRâ~šÞâJñ_‡ðÈÕ" },
	{ "rong",L"ÈØÈÛ‹†šÕéFŽVÈÜ˜xÈÚÁsÈÖ“rÈÝžq¿^ó“Ï”‹’È×Ñ’ˆcægÈÞ‚æÈÙt˜sëÀ‹æÝP½qéÅ·ZñŒÆŽ•íáõ¬ŒòîÎ hÈßŸVáÉ‚ÔÊŒ]“mš¿" },
	{ "rou",L"÷·óÎj…œ¬yŒ`¶bõåÝŠÄ\ù’ök‹Y»€ÈâôÛÈàòkœn˜QÈ|åˆíq­~ÈáŸ§" },
	{ "ru",L"á}ß”JÈìñàÈëÃNï¨†äãœŽš•ãÈè–ôå¦¿dò¬ÈãÞ¸Ýê‚¢Ê‡ÈäÈåÉSû«Aàr ^œx…ÊÈéîž‹‡ÑMø›¹TçÈ’äáÈçøn÷pàéÈæÈêè`Äžä²‹çÀ]" },
	{ "ruan",L"‹\ëÃÄQÝ‰ÈíÈîÜ›Þ ^µO™“‰¼¾¬}Îp“ÉÂX‚¢­wˆë´M" },
	{ "rui",L"ÊtèÄÎT»ÜÇÀB®cÈïÈðä„âcî£¾q…±—MÞ¨ÌGÌH™G›I‰ÇÈñÆò¸ŠñäJ" },
	{ "run",L"étéc˜ôÈó™Èò" },
	{ "ruo",L"Ÿxóè…ª‹SÙ¼àeö}“ÉÈõö”úUºOÈôœc’Ú’µ—í kÉm" },
	{ "sa",L"è•ïSØíìƒž¢çoìªæp›‡È÷ì‘êýëÛË_ÔQØ¦ÈöëM”cÜaâlÊ”ñ`™¨¥Èø“—" },
	{ "sai",L"†ðàº›î|ÙšËÈúÈûÈùÈüË¼šºàç‡Tƒwöw“H" },
	{ "san",L"Ž¥éd‚ðð€ôLáê‚ãçD ÑôÖ—Ø™V¼W¼B…£ãßqâÌ¿™çoÈþ‚^¼VÉ¡šÉÈýÉ¢ë§¼RšÐ…x¥" },
	{ "sang",L"†ÊÑ˜É¥òªÉ¤î‹–øÉ£æríßÞú˜š" },
	{ "sao",L"öþÉ©ÏAðþ¿„ëýš×É§É¦ö…çÒò}²„ÀR‘¨ïbçØòX’ßà“ýèA¿‰œÐÉ¨ó÷fÉÒÜ£" },
	{ "se",L"æšoÀNãG®æaï¤ZÞQšmÉªÉ¬ë–ÜäC››“ö’‘‘­·wêSÉ«Ö ­·†»¯™†Ý¬Xð£ÈûžiØÄœàïo­içm" },
	{ "sen",L"ºdÉ­—Ø˜¦ÒI" },
	{ "seng",L"É®ôO" },
	{ "sha",L"ôÄƒƒ’­š¢ï¡É¯õçéŒì¦—EŽ¨àÄ»}ðð˜×³À\ßþæ|ö®¹€†ÃëÇÉ³¼†ô‹“—„xÏÃÈSé„É¼˜fÉ´É·É±É°†~õÁ †—öè‚êýÉ¶É²BÁœªQÉµÊe" },
	{ "shai",L"•ñÉ«õ§ºkÖLÀ\ºY»iá‡É¸É¹º”ƒ" },
	{ "shan",L"¸ž¯ZÉÃ’ïäúž¨Áƒëþ×iáŸÉ¿“ªGÉ½’´‰æÓ„šØßþŠê„ÉÇÛïÒv™c–Å·_‡AÉÈ”»Ü‘ìøðƒ¿ƒ{˜èãˆŽE±˜Ãˆ÷g—Öæ©™ÒÓ@ðÞ²ô“Û¸–É¾‰ŽÜÏ¶UÕ¤µ§ŸÄóµÉÀÒIßò~ÉÁ•Ï€¿„Ó˜ÉÉõÇô®õŠŽ»éWµ¥þˆ†ÎîÌ÷WÙ ÷XÉ»Ú]‚ÞæóîtÉÅø@ÉºƒR“½Ê`É¼ç—´ŠŸšÚ¨Û· „hÉÂ÷­Š™ÖbÉÆî¿˜ŸSéXÁÀuÉÄžè”vþ…" },
	{ "shang",L"ç´ÉÊéäìØš‘ÉÍAÚJôlœ«ÌÀÉÑCˆsÕõüÉÐÉËÙp¾yÏDÉÏv¶@èlÖ…g‘^ÓxÉÌÐLÛð‚ûÊKˆÃ‰j‘ûÉÎŒ¬ì " },
	{ "shao",L"ÇzÉÜŸýÈpÜæÔ½BÉÚŸ†ïYÉÕÉÙÈVímõ}”ïäûÛ¿ÉÖÐŒòÙó™–¶¾Kô¹É×è¼ÉÓŠ¾…pÇÊÝi±ÉÔÉØÉÛóâÕÙÊ–«x„ÉÒ" },
	{ "she",L"Ï‡Ýf’¡ÞéÉÞÉèòM®ŒÙÜ‘ØÊ°ž—êA™œhÉä‘bíH’wÉá“ºì¨ÉçÉàÊJÉãê^’Î™Ýî´ÍF›õØÇÅhÉâÙhÔOÕÛísäÜâ¦ŠL÷êÙdÍ…”zÉæ´’ÉßÉÝ…‡Éå" },
	{ "shei",L"ÕlË­" },
	{ "shen",L"Éð•Öò×“äÉ–¸ÔBÉôô…¢Ô–•YÉêé©¾D˜¦BÃŒäv»r²sÉ†…£‹ðvù_Ÿö·ŒÉóŒqÓ\šáÈ²ÎµŠÉï×}ÉíêÕ”’JÚÅëÏ÷“¼R±mõ˜ãh’bÉöß•ÉñÉëŒ›Øñ‘õÉîÉò»põŠŽ»Š·ïò¼ÁKöŸ«|”žÊ²Í–Éø·ÄIôÖÉìÉõƒÂÊQ•Ü®e‚L¯”±s×Ÿ¯}Ñ[² ²_‹Ž—ªÔYœVÁAöYÉ÷zžcËMÉé…¤ˆÞ®`Ý·ßÓŒæÚ·ÝØîTŠ" },
	{ "sheng",L"êj õÊ¡„Ù…Ö˜|›ˆÉùš}ÉüÆÊ£Â}š ù|œ¤Éþ‰˜‚¯Ù‹ãHäÅ”Îå•ê…¿IÉýê’„œƒÉûÊo¬]Š¿íò«{×WÀKÂ•–™ÙKÉúÊ¢÷jÊ¤•N¸i‘™•úü›Ÿ„Æ•…Ê¥³Ë\óÏáÓêÉ" },
	{ "shi",L"‹ÒŒp‹qŽŸª{•EÊ­Ñ öõÒ|Ê©âÊÔŒÊ½ÖÅ«ïz•rÊÁÕžýaÛõÊ·Ê¬ÊË†Fœ¢±c¯aÊ¨•gÛÊ°ÎtÀ[ÊÇõZŒÆŠ]¶ƒÊ±ï†á‡â‹‘÷dÔ‡ÅkÊÑÓl…„ìÂãAœÒï—ÊÎö|ÛJ…áøOœ›áyªL»iüœÊ¦ãJ³×Ò•ñ‚Îgöå‰PC˜tÊÃ±xŒgóÂŸ³ÊÉÊÏ½J^–É–§Ñ|ÊÄ¹EÊÐõ¹ðOâ»öˆ‡ð™yÖ³ñÊÓâPÊ¹¸bþ˜±Ònó§Š¸ÝªsîæÔŠ]ðSËÆêÛ¬‹ÝéÕœœÛ~ÚÖ…«Ê¯åœóßÉNÊºÊ¿ßrÖÊ³–òƒ½ÉÊÆßŸû\ßYÊÀÊÌö‰âž¹•ô“JöXÉP›¸úPß}†çÊ¶ªHÊ¾õ§ÊÂéøÙBÊÒÊ§ÊÊF±iÊÊ¼ÊÈ„ÝžøœáÊ²ßfßmÈžÊ®×Rƒàß±ÊÅŒjáŒø[”—üô˜ÞyºìêŒ‡uÖuãvÊµµuäK…Ú…bÌÊ´Êªñ\ÉÝYÊ»ÊÍÊ¸Ê«Ðê" },
	{ "shou",L"ˆ–¯lÊÚ¾R«FæÊÖÊÜô¼ÄfÊÝ”™ÊØÊÛÊÙÊ×Êì‰ÞÊÞ…§á~‰Ûç·›ìá÷ÊÕ" },
	{ "shu",L"ñâÊåÚHƒ©ÊôÊüÊßÒl­qÛSÊã‚T¼^àg–µùŽƒ¬ŸêxË\ÇO¸w–XÐO‚måf¯EÊë½RÊáÁ›¼‚ÑVÊ÷ç£˜Ð»PùeØQçTì¯Ø­Êâ–€ïøÉDÊêÝÄÊûÊðÊþŒÙÊö”µ•¤’¼ÊøŠì¶•†CšÑ”d˜äòÜ“Œ«©úžÊõÉ[Êñˆ™]ÊæÞóóŸYÊé°P’¿âàÊìÏŽÛÓ½ˆ÷t÷nã_Û\XÒeËŸÊçæ­ÊèäøÐWÊòÊúÊàêœƒÊÝ”ÊýŒ¥Ù¿Ì ëò‚JÊíË¡ë¨ÓáÖ‘Êï•øž‚Êäãðõ_ÊóŽõšÌ’æ‚‚ÊùÊîÐgŒF" },
	{ "shua",L"à§Ë£ßxÑ¡Ë¢ÕX" },
	{ "shuai",L"Ë¤Ë¦Ë¥Ž›ó°Ë§ÂÊÀŠ…i¿\" },
	{ "shuan",L"Ë©Ë¨éVÄYãÅ˜¤äÌ" },
	{ "shuang",L"wž“™ÜÆCç`‚öæ×ž{ëpût‘SûUú{¿Yò‚˜¾Ë«ËªóLãñË¬‰u‹þµdóZœö" },
	{ "shui",L"µˆ¶šìË®ŠÜË¯Ëµéj›çÕh’¨ÑcË°Ë­›äÕfÃŸ’ÉÕlŽœ" },
	{ "shun",L"ÊŠË³Ë²Ë±ÝôB²p˜ù²iË´±†í˜˜J" },
	{ "shuo",L"´Tª“Ë¸Ë¶ÊýËŽÝôÕhîåèpË·Õf²æl dšF”µËµ qéÃ¹››«åùÞ÷†d" },
	{ "si",L"‹w‡z´fÀŸÌï›qÌŒæ¦ËÆñê¶Dâ–ãôƒÞ ïÈ—t‚Ð¾ŒÎ‡ï\Î’æá–yÙîºÒ–âL˜{ØËÙ¹²Þ›…Ø|Ê³ìëôé›—ËÃ¶Lð¸Ëº—öñ†ý½zÏa¶TÛÌËÈâ‘çÁÁQãáË¹¸r µòlòIäùàÃB‚hËÂ–ŸÊœòÏÖpçrŒKãƒž[Ë¿ït…‹PŸùÎEË¼ï~Ë›ålæJ–ÆÏz„@ï•Ë¾úƒË½Ê‘ãjËÀËÇËÅóÓälúfËÁßÐËÄýDË»äFŠÙ‚Æ" },
	{ "song",L"áÔÚ¡˜BÝ¿ËÉ’¿Î@Y‘¡ËÊ‘mæ–·—ŒËË–…äÁ™€ÕbÔAØþdËÍëþeËÏ³—ížñµ“¡‚öðmËÎ‘ZÂ–ŽôŠ»Ìtñžè—sËÌßâì‚‘ËÐã¤áÂó " },
	{ "sou",L"òp…®æC“¡âÈàÕƒðÞ´ËÑ™¸”\ànäÑî¤‚ÏËÓ—¯ðtágÉrÛÅªvïËšFì¬ïb“–ËÒæ}ËÔòô’ÈŽù¯˜Ë’ï`ÉL»Pà²" },
	{ "su",L"¿sËÝ˜ÂóXË‚Û‘­X˜j«T·dä_«Žš÷T®d¸@ËÙ‚ÑËÞóùˆ¼ËÖË×à¼ãº™ÅÌK‡Õûh…r˜þ—­ËØßi‰OËàËßœßÇxÄhÃCä³ú‰ÔVËÕ›ƒ—Vö¢ÝøšƒðMÚx¿iõ‡˜É‘ˆÖq»´cßpËÜËõËÚ‹•’ò“ÙíÌVåËÛ›«ÚÕöÕ" },
	{ "suan",L"¯i¹gËáËã¸Œâ¡µ{ßxÑ¡Ëâ" },
	{ "sui",L"ŠÌœñëmËç·ußU½—Ãœì[Ëîçw™p·[êyÀŠËäËæ­…Ò`¿\Ç]ëS‰å†a†÷Ÿ«ˆ¼ÀšËŸÕ¿…å¡ç›Ëë’µ³ZíõÕrÜ‚‹Ó×\ËêólÙw”øÀZËíÄòÆVËìÝ´‘Îí}ÄŽ²Bçiåä‚‹›ÔËå¿“­j¶XËéËèìÝºwìšÚÇÈššqžvî¡" },
	{ "sun",L"ËVÉpËïé¾ŒOúZº‹¹Sªsáøâ¸ïŠËñÊ˜“q†Ð–Ëð“˜æ{“pÝ¥–Õ¹öÀ˜ƒ" },
	{ "suo",L"Ñ–ËùÏ×’­¿sËôæaíüÚtÉ¯Çjóš«I»‚é•­ËöË÷õ€ßïèø ÞËøôÈæ•­Fæ¶Ëóæ\ËõàÂºwêýºz†îàÊ‹‘“™œÅËò¬æiô‹æ" },
	{ "ta",L"“Ò‡Åêê[ìŸËý«H™\ßQíOœÍêY›øê`¶N÷£“‚ædãËÕw‡–ÍØ…úí³ ­Ì¤äâæ]šÏäðËûËþåÝÌ£‰‡Ëüí^þjÌ¢Ëú˜dõ]ãBé½õÁ‚èÌ¡ÇEêFîèÜcß“éß_„ßeµkàªåJÑöªH×n" },
	{ "tai",L"kïUÌ§Ì¦™…Ü–”E–Ÿ‰ûöØ”ÁÌªîÑÅvÇ œÌ»FâÌ­Û¢Ì¨Å_‹ê M…õžåìÆƒˆÌ«ñ~õÌŠU‡ò¹xß¾Ì¥Ì¬‘Bƒè•@«}ææõTÞ·ŒLëÄÌ©çˆrö" },
	{ "tan",L"Ì±ÀW‚„‰¯ØÌ»ãgÅjðZÀ—‰ Ì¹ˆÅá]Ì²å£êæ¯a‡d‰›ÌµÅl‚èÌ¯µ¯°cƒNîãáaØ”Z•Òç†Šò‘˜ïÄÒfÌ³”‚ñûáv´Ì¸ê¼šU‘…Û°‘Ÿ“ÚÌ¿‡@×Zœž­fž©—Ì·Ì¾îtÙyÈIÌ½ú‚Ë“×T•Æ†ú‡cÌ¼Z ˜W¾gÌ¶†®håU‰ÂÀìþ™AÌ°Ì®Õ„Ì´Ìºïâ" },
	{ "tang",L"ºLÌÈñíôÊ CÌÇíUœ«èKç|ó¥‚«‘Üæ†ü‘´g‡RàûðyÌÃêOçMgÛ˜üßT“­•òf²˜ï¦àoéEÌÂæh”†ðh†°õ±ÌÄÊŽ¼CÉ¶KéÌº‚è’ÌÊÌËïÛÌÀÌÆÌÅÌÌŸ¶äç¹ÚZ˜yêWÄgúSÌoüh™éó«è©¼QËTÙÎé‹ÌÉâ¼ƒ¯ ‡ëGÉy„¨ÎvˆnÛ}ðnÌÁ‚Ú" },
	{ "tao",L"ÌÒÎIþ¬•ÞŽµÌÖ¿l†GÌÔÌÏÌÕ‹—ÖzíwÑi¾Iß¶Ý‰ú còP“†ßûðu„üÌÓ—˜…ÌÑñŠ¾TÌ×ÌÎìŠï‘ÌÐ™„À‡íN|èºýÓ‘ÌÍá[—ƒì’ä¬åcØ»¿_½dìâµŽ÷Ò»Iä•Ô|›ì" },
	{ "te",L"ÍfìýÙJÏcß¯í«Ø–ï«ÎŸäˆÃŽÌØ Ã" },
	{ "teng",L"Ä†Ö`ÌÛü’¯\ëøóIÌ„þ—ƒ£öŒbŽ¸ÌÜÎŸƒ\¿gÌÚß‚»LÌÙ»Tòvñ" },
	{ "ti",L"´f–õÌëŠ¸‹qÁHã©šÌÝ“WÌŒÆl¨ÌÞÌéÚ„˜Nó›‹X½úeç°ÜnóƒÐ}”`äRÜè¶”ßXóe‚mów›¢ÏsÌàÞ…Ìãö[õ®òf†ÙóÌæ¶_ùY‚¨ÖBÌä‡¢wç¾Ó¾ŸŠ­ƒÙÃåÑþÊø˜ŒÑÖpÞ‡Üƒ…†zõ{ËSóžç‘«Ÿ¬vÊƒÚŒ»G¶Aù•ú‚ŠDåÜSÛ‡î}·aÌßœvñÓ’«Ñ|ÔgÎyÌçúfÌâÌèµðÃÉ‘ø÷–ÞÐõkšYÌáÑ{ßPŒÏ’óÌåÌê´Yù—" },
	{ "tian",L"ŠÇüV“¸Kå`Ìð¯t²V¬™úcâšã”œï›ìtŒÄƒÌáLî­k¤ÃbãÙÌò‚’×Ùq´kÌó®\Œ…†Š‰\±]ÞÝŠõÌñ±™Ìì•‹œLîäÓ`úlÌîìpÌï›pÅqî±µèéå²_îŒêDèÓC¾gìj´[¬_åU®ƒÌíãÃÈJ" },
	{ "tiao",L"î\ýfÉ‰ì›÷ØõÆKìöµxòè†G¸IÜæÚqÉŠö¶—lÌ÷ÌõÙ¬ÈVÌø˜Ôföœ¼gµ÷äpGÏC–Iñ»ÌöÕA”ÓôÐŒýæxÕ{Ò›•qÃxã“ŽçÌô½róÔÂwŒi‹àÉ‚”þÅ—ƒ©öæ" },
	{ "tie",L"Ìù…ãâŸãŽÛ@÷Ñƒcø‡èFòÙNÂzï”ÝÆÂç“ÍuÌúä~ÌûG" },
	{ "ting",L"ÂŠèèŠÇ…ˆÃ‰ìŸP—þ¹j¬EÜðîcÖF›àì˜ŽØæÃöªïFòÑéƒî®‚KÂ[ÆJÎbß‹‚DÌüÍ¤Í¡Í¥ŽßaÍ§ÕPÌþ«žîú…ÝãâÌýüžÍ¢˜w½–ŸNäbdˆNµÂŸÂ Í¦Â——HœsÍ£" },
	{ "tong",L"á¼½ŠÙ×†LãnÙÚäüÍ«Í°Ù¡›Ïã~¯]ŸüÍ¨žú·rÚU ‚Í®–SÔ˜Žä‘q¶‚•zÛíÄ€ÐhÍ¬dÍ¯Øçõj»Í²ÉŒšÔï ª‘½yÍ±Í©‡ìâú³‹¶²Ÿ×àÌM˜¿Íª½pãPÍ­÷‹Í³íÅ„ç Õ¹cžç¶±ÍU‚£Í´‘Q„¨™HªIÜí±•Ó" },
	{ "tou",L"ˆÇ÷»Í¸äWî^Ùï¾–†V‚ÊÑˆ‹ÍµÍ¶æBÍ·‹U½‘" },
	{ "tu",L"îÊLÍ»Ü¢Œ_ùIå„ÛT¯…ú“ƒ·„ÍÃxùrœ£ñG¬ŸŽêõ©ù“\äŒ’ØˆMÍÁÚ¢Ä‡í›ÞÄ]Í½Ý±âŠùW›BýCÍÂÍÀ†lÍ¼¯fÍ¿È‹Oá“ŸúhÚgˆE¹\Þƒ¶dƒòÍ¾âQÇ¤’¼ˆàÝË¶•ÍºÍ¹òB—^‰TÉ\ˆD†ž" },
	{ "tuan",L"‡âˆC®™`ù‡É”Ñƒ™ˆ÷HÞÒ¼aœ¨åèØ‡ˆFÍÅ˜¤ªl“»ºiú™„Œ‰’ÍÄÑ‰„–‰tæ˜î¶‹§Ÿ™úo‘_" },
	{ "tui",L"‰‘Ã•É—˜úß¯ÛƒÍÊÍÉÂv×‚ëPîj·~Œ¾‚Q‚MŠÑôsÎ†ƒUìÕÛóhîkòoË”ÌLÍËòDÍÆŸlÍÇÍÈînÍ‘" },
	{ "tun",L"â½Ü”ë˜®™ÆXÍÊ÷ü`…×•HÍÍÎP›â÷ƒÍÌ–Nï‚ÍÎØZÄ†ëàÙÛ¶ÚêÕÄ™ŽÝŸl†”ˆd‡pô" },
	{ "tuo",L"š¼õDÃ“Û|’¨Ø±’ÉÈ[‚Mô…ÍÐ—‡îè öö¾ÍÑ½FšAãõÞ~ÍØ÷WÙ¢ÍÖør“ã‹µâÕóCÇhãB…ïñXÓšÚ—˜’ÍÕ’„óêÔqéÒÜ€×™ò™ï€üƒñ„ü˜ÍÙñWözñ…ÍÒíÈõ¢ëð˜ÍÓ±‹sèÞ–lÞÆÇõÉêu»X ­ˆ÷šÍÍÔÍ×—ø›ñÍÏêeãûÌE›kèØÐ†ñjÐ›šú¶æ³aÛç’L´P™EùK" },
	{ "wa",L"ÍàÍÞœÎÍß„¾®|j…÷·“°¼æ´ÍÚÍÛ†åÂvëðŠ¹ícßœÍÝ†ì³[Øô‹zÖœíiÄeÒmŒÜ“‰ü|¸Dí€†„”…º†œ·˜ÍÜ·Šì…Ž’" },
	{ "wai",L"žxáË†JØ²Íáî“†·þZÍâ¸" },
	{ "wan",L"ëäÍèä‘æýÍí‰îÍñÍðÍîà„¾Oçº•ŠÍéŸêPÍóÃäŠ€Šþçþ±DÍäÂû›ð†n²o±›ÏÚ@ÍçÈf˜´’Âó[Ãévó\Ûl­ÍëÙ–—µž³ÍæåsÇ|¾Uå†Ø™—iÝnËHîµ’eÍêØžÍì‰GØà‰Ïóî•–ñš÷„\Ý¸äjŽ¦’Ìòê¬T‚{…dô’…eÍòÈxÝkÝÒÏT•ˆÂDŒñÜ¹ó]ÍWˆ¾êK¼wÍåä[îBæ~‰íÍãÍï" },
	{ "wang",L"’[À Íõ±ZÍ‡ÍüÍ÷Œ¶ƒÇž_Ýyã¯Á@ÍôÍùÎ\ÈDŒµÍú¬]—Ÿ¾WØèÞ‚ÕsÍ^÷Í•™ÇwÍýŒ·ÍøûÍûÞŽºéþÍö–R¸Œ²´“©" },
	{ "wei",L"Î¸÷˜Î£Ó}Ÿ˜Î¯ŽUÈ”õnÎ¶È–ƒ¤Ÿ£Ç‹Ä^îQ¾•Îµðô˜LžéàŒ¾Sï]ÎªÁWðŠàíÎO†ªcÝÚeä¶é°IÆYÎÓ„”Íãí‹WåMÎoæ¸Õ†³}ŸÝ†Âó[“Öðj°LÎ¼­Ž ÒÜZít³u‹n†Òå…ó\Æ„Î¬ÒE‘Î¢—Ü…yìS“ã—ÛÄŽ›¾ÓWÎ½Î¦è¸Î¿ÓAêžâ¬ŒÎ»Î«â«õdä¢Ëeçâ‚ÎÙËìGÞE VœwÌå—Î°œ‘‡úìÐ‘£ãÇáË´SÎ­ÇUê¦ÏGÔ•Î±÷Îƒ^Ü^Ó‚Îº‰ÃßzÎ²ÉJáW¾“Î·özÎ¹ög›W ‘Ú~ÚóÛc—|´oÞ±á¡×ˆžw•¥ÚñöÛžHÒÅžùÎ¡Û×@¯_×~Î´É–Î®Íþ´j¬ÊlË—Î¥Î§›”Ž®¬|“GÎ¾­M‹yÎ³ílí|Î©ÛbÚÃ¬^áÍÎ¤ß`žS_ó]œ¿“fŒËð]çAôºœÕõKèÖ^åÔ‚¥ÐoÎVàø à™ÌvÎ¨iíf‰Šì¿ÎÀÎk…°Ðl’Ë½öh" },
	{ "wen",L"Ê•…Ðö©ÎÄãëøjÀˆØØšzÃ‚ÎÉžÉ“hÆ[‰eËœ÷—œbÎÊëÃWÍPÙï·€Šp¨øs½ƒéÅgñbÎÁÌNñmü•Ý˜ÎÇÊŸ¿A†–ô•˜XÃÎÈ¬ÑŽ“‹ÎêZøY¿Z¼y—S•jæ’’^ìíy¯‡éš˜v…ØÎÂšØnÎÅî‚ÏRé”…ÝÔÌãÓö€ÎÃÝœ·gÞdÇ|Â„ö“è·œØÎÆ«œ" },
	{ "weng",L"ûl‰R²\Ç•²úO„ØŠTÀšœåæf”wÂýNÎËÎŠ®YÎÍÞ³ÎÌÝî" },
	{ "wo",L"íÒ›óÎÒÎÔ—çà¸“ëä×œuëo¥ÎÓýŠ–†ÅŸ‡—Îá¢ÛlÎÎÛb’Üüªi­xÎÕë¿ã’ÓËhñNý}Ÿs†É^ÄŸ‚¬ÎÏÅPŠðÝ«ÄOÎÐÎÑ›ðˆå¸CÙÁ‹_Èn’Úö»²YÎÖŠñ" },
	{ "wu",L"âèÎãv†•‚—Õ`Ïwžõâä““ßíÎñöƒàw…Ç„ÕðÍ‡f»|ßAŸÊì¶“h…Ò’HÎ×âEŸoåüän„–fú~×OÎèìFíÎéº«bÎîúFÆ•Îï’GÜÌå»èžµÎåÎ‘“õˆÎäë‰ä´Ž–gœrøŒÎçÎÝšTƒÎì‰]ýIÎÚÎÜÎë›^òÚšœ×˜î”–Ç`›´ðÄ’æÄìW}ŠÓ‚W›@ýrÎÎðŸ½ýH“ŒäÕ_ìÉÕGŠÃÎê²y¶³Jäoò\ù^åÃì}ƒÇêõÎóTÎØÎæÜRÎÙ¹™†èÎíÚùÎÛê‚ÄŠè»ûc÷ùëœÛØ°ÎòðíŠV’N¶ñŽÄ®W·—¬@ÎÞÎà« Ø£ŒíÊÎá¸P´IÚãæu›A‹³ëFâÐæð•J‰­NÚ›žŠÕÌFÎâØõµŸÎßÍöàNùMòú¬öÈ" },
	{ "xi",L"„DŽQ½”ÔD¸OÚTÎùÖìùÏ¨Ž`Ï´Ï²ŸÁŸè–y‡q“©¿ŽãÒÏª†ÕãcÌŸ•ÊÜhÏ¶òá˜~ë^àEšãìUÎôØlË@‡½ÁpØ‰ùTŽdÒuæˆÃ|ŒÚjÆœë•KÎýõèßñÐÏ¡±_ Þè•¦Ï«‘‚‰Iá‘ïÀ°ÞÉÚv‡ÖÝßÏk´mô¸ðFª“ŸXåïŸm•‘ÃZÓ‚‚ÝÁ—æÒìäÕOÏ¬Î‰ÊDÆÜÒ ðœÈ}âR­Œì¨´—Ï£ŒÈ­tØgÃ[ØGêêÏ±óNÛ­ÏµÇmò„™ú÷žàSÉtÏ­À{Ï¯Û’…ÀŒjÀ¯Ï¤ÎEÀðqÙâÁ•â¾ÎûÀG_ôªëKìûÚVš]ŸyýA÷@Ó}ÐPÖLÎú‘òÎüÇbŸ_Îõô]¼Y‰€â|ä»‹ÄOÎþÎöó£‘ñðª‚S’VÏ©ó¬…c¼š•„×@²qdÓ„ú —ÌÏ¸Ï·´Fïôâ˜›Ï³ªLÃ~åa¬NÏ¢äÀ÷ûÏ°÷^ÉYŸùÎ€âl xŒÁÏ§ÖlYïeí’Q‚`Ê“ì¤ñÓòwŸ¼À…ôÑÓBÏ®…äðOÜç˜éá@õµÀMš@Ÿç÷ŒÊÚiÙÒ¾küŸë¿u¹Ï¥áãêSÝû‰¸– †Œü_ç^†AØHÒ‚ñ¶éØãbãŠÎ÷±–ŠÖ Ì´ŽS…kÝ¾¯Œ·G…wÅbêØÐaìIÎøÚÀôËçô—áÏ¦àq—N”Úôâéø¿]ßëvè„Û§èÉj‘ƒ O²—âMŸ›…[¿" },
	{ "xia",L"æ_ÏÄ”¯‰ì¸—÷ï{Ï¾ßÈ’Üç]Ï¿å’ªM¹d•gÏ»ÏÅÏÁžþ›ÑŸžÙÎ˜ÏÀÝ íÌÁŽÆS…­èÔ“Šô öyÎràAòhÅrË´WšB³ˆåÚáò‚ÒÏºÏÃBÊ›ˆ®Ï½óÁÅ{æ‘³´lµ„ÖlÚYïPêƒ‚bÏ¹Ï¼¿[‡˜ØB«”ê˜è¦BÏÂ¯KÕ’»£²LúT" },
	{ "xian",L"¼`óÚÐjÍpéf†é½ÝËÏÌ°Gí„üí`ÏÇÌ_²v‹¸Ï´“ÍŽÒÏÏ½žë¯¿”gµ »˜ØRå‚µUƒmÇ{õÐñMçoƒMýEªž«ˆÀo«NÝü†m›×ƒgðW‘œã•±]ú’±•ûy·SïÄ‹MÒŠ†Zées½ÙÌ`Ï³ÏØÃjÃ|Ã~Ý²†¥ÜŽÁáýÏßÖ›ŠÒ°B¬FÜ]Í€ÏÍ¾Q¾€ÒDèËWÕtîy˜óá_ÏÒÏÓåU±hÏÝÏÜ‹ÍÁ{ÏÖ–}”sÛŸõrŠÞ³wï@¿h‚]ž¶Ï×ú‘Š½‘—õ£Ì\‘`«IŒ¯ã”žnªA½måv…ûÚDœ¶…î½LÏÔÅ@“{ÏÈÒv’¦ÖP“ÈÄdÚ`á­Áwò¹ðïðÂªÏÞ‰·Æx™ÌŽMÏÛÏÆëUüGÏÎÏÙÀ‰™Š«ÏËŒÝììÏÑúšãŠšÀåßÅ`™÷ìÞö±¼û¿„ôÌŠ·ÏÕžóÙþ‘¾Õ^æµÙt²{‹¹ƒn•÷€ÏÐÍ˜í†˜åDŠhúNÜÈû’®QÏÉ½Œ°ÜŒ‡JÏÊÞºÅOˆŸÏÚ×]¶iÀwä}¹‘õÑ‹ü" },
	{ "xiang",L"é•½|Ïâ•}Ïçí‘ç½Ñóã}½µàmç}Ïåø—ÝÙ–Ùãà_ÏáóJè‚Š¢ÄÁfŽûÏãÛKÏñÏïÞ†Ïó¾|«“àlÏìßðA÷ÏâÔð“æø÷zàxÏòÏéÏæËG‹Ïðô\À‘­ó­Ô”õœÏ†ÏêÌZ„ðÜ¼†“ÍJ·E™ÖÏî‰Ði­˜‡»Ïä÷Pƒ¨ÝâÃñ•ÚÏëí—û‘é{ÒVð‹ÏèÈeößé—ÏàÏí„âõaÀv" },
	{ "xiao",L"·nÏûóïÏvÇzš¥›©š®·›ÞB“ßïYÏSŸ^³‡”Ã½‹Õ[†Û“`ž¼†DÐ¥Ô‰Á›û^áÅúr‡[ª’ª”ÖjšRžñË@ŠëÐDŽéø“àUkóuÌŸêÏöòÙÏ]•š‘‹÷ÌÏ÷ßØ‡EËrÖy‡Æújø{±Ì‡èÕÐ¦Ð£šYŸÀ–žtx™Ïº}ónš^‡^Ïü‡VäNÈp°~Ã‘‡Ì¹qºSÐ¡°†äìÏúÖ—”¬ž½ÕqÏý˜þ‚jò”Ð¢æç•Ô‚åÄ…Ê’ÏþŒnÏõ _Ð§”Â›ßóãåÐÍ“Ï…ëÏù‚PèÉºÏø¯hŸò¯eÏôÔF‡CªVøŸÐ¤„¿—n×Dç¯ÛX" },
	{ "xie",L"ÐªƒDæE”Xí“ýklÐ±•»eýaÐ¹µmÐ¨ôkÒp¬€ËZÊÐ³‹rç¥À‹¼I”ýöÙõ@ÖCâ³ÙÉÄnÀTß¢í…–¤ŠGÃ|Ã~ŒÚœëÑ€Öxï¾Š„µÓi½âˆ•ŒÁÐ¸Œ@íP‚´”yŸcÒ¶éÇý^Ð®‰êÐ¬é¿ÀiÑªÐ°×Ð´ÆõªnÐ~å¬…lžaŸLÐ©Œ‘ÒCÐ¼½uÐ»Ò³¼œÐ·ÎdÕ™íC½’‚ÄÃ{ýš Xžá’¶ÇäÍŒÔžÂ¾™½X†àèHÏÐ¯ŒÈÐ«Ð­œœ‡ƒ½eåâõóâÝ“ûƒª¶cÙôŒÑ“a“y yõqË†žàŸ»ÜaÐºò¡‰fÐµÞ¯ÛÄÁ–Ï’´cø…fÄîRƒæÎqýKÐ¶ÛÆŠÀ’çÓŽO›ªÐ²" },
	{ "xin",L"âdQþ€ÔDÃ’ì§ŒJÐÅÐÂÜ°–“ñQÒWžÔ‚r‡Œ‹×Ø¶ÐÃ–‚êcá…ß”îˆ±^Š|Ð½ŽßµUÐÄôgÐ¾¹ä\ÜŒöÎ¿Ð¿ñ^Ÿ{‘€ïâÐÁÝ·ê¿ÅgÔMç†ÐÆÐÀ" },
	{ "xing",L"ŸÉ õÊ¡â]óU‹”ÐÈÐËòHÐÏ¹“ÅBè™²M‰D‚†ß©ÍÖ_ät•Û›ëˆlÐÊõSÐÐÐÎè—ÐÓ¬wÚêÜþðh¹žÐÇ‹ñÓwíÊÐÍÅdœîãoã¬ÜôÐÑ›™ŠüˆžÐÔŸ“Óqé°‹ŽyÇnÐÒÐÉö]ê€ÐÕâ¼ŠÈ¾mÐÌ" },
	{ "xiong",L"ÐØÐÛ×œÐÙÔKÐ×ƒ´Ÿ‚ÐÚÔwúÔžÜº›°×›ÚUÐÜÃr„öÙ‚Ÿ‡†Mr‰éÐÖ" },
	{ "xiu",L"Ðå­Pó…ÎÉŠÃ‘Ðß½‘ÅW÷ÛßÝ€«‹ÆvLËÞ‡›äP÷GžòÐáÐÝ¬Lõ÷ð}žñïqæ™äå˜¼ð¼³ôÀCçn‚cÑ„Ðã½œúÐÞã–Ðäø ÐàÅ^ýM¼NÐâÆ’á¶âÓÃƒâÊÑ…çVæT" },
	{ "xu",L"š_…•vÐî¯LÐõží›U‚Tš~¿H›TÐ÷Öž¾{ñãÔSäÓ“TŠ˜ã„Ž­„Ô«Ô[Â{±S—ì‚»É’Æ^Ðô”¢ìãÖ[Vœ•çïGíšÕš™øªšHÌôÚÐø•B¾–±r‰ÙÊŒ‘AÞ£É[ä°• Ðæíì‹€ôzÓ’ŸTÐéÎd·PôqÂ…ã_•d±Ì“Å²xÐìsÐïÐðË…ƒÛÍ‚…räª‡bÚ¼Û×±Ní¹†ÄÐçÐíÐó²W÷ròš[ˆ¦ÓõÐë× Ö~ÐöÐò•ýšAËvÀmœä‹ÁÔ‚õ¯‡uôPÙ[Ðè«—…é· ’îÛÃ”›Ðñíœè`èò·V¾wÑSœMà†ÐêòÀ]" },
	{ "xuan",L"ÐþÌBäöÈkÈ¯½kÐýÝæÊžñ¿’™e“Eé¸ßxìÅ•ÃämòCª™ÞïìÓ²Ñ£«Rè¯•RÂA‘¤…éIÐfÑ¤×XÑ¢Ÿ@ÜŽæM²U•]ß€¬KðçÏŸœÊR°_äÖ‰H±†ÙØÐúÍ›¿hèGÐûÞFÑ¡ÚK…ºêÑÌTæ›Šˆ•œ½LÐùãCÚÎÐž¶PÐü‹Öîç‹lÉ{—]×‰éÂQ­vË‹ŸÖXìœÕíÛ¿•tãùÎh¬I‘ÒÍ•˜CªBïX¹ŽR†IÐ«t¬uïà" },
	{ "xue",L"–ùõ½ˆyÑ¦›‰ú›û`Æ‹l…ÉÚÊÑ§ÚpÞm²xÓ{ üžyŒWøŠNŽGíYÖoÉHëzV¯TÑ¨Í‰®œéí| KÑ©”ÄÞjÐÏ÷í´›‡÷¨Ñ¥àåŒú÷LÑª¯N" },
	{ "xun",L"ˆ_Ïr oÑ²¶½kñÑ¬¡­RèRŒOÑ¸˜ßÓ–Ÿ[Ó¼rªF„ëÑ´ô‚Å…_ÈÑ«ä­ßdŠQõ¸Àc‡x¿£Ôƒâ´ÏyñZš¦‰_÷\žF–hÑµe–Õ‡ ¾Þ™‡eÑ­à‰ñ¿²†„×Ë`Ñ¯ÙãÑ¶Ñ³Ñ®êÖÑ¤ @ w”Þ¦‰¶Ý¡Û¨ÒWŽ…Û÷öàÑ±žµŸñ»çÞ¹âþÄ“MîšŒ¤Ùbä±Ü÷á¾ÌQÑ°Ñ· `š½„ì‰Ëù—D®p«‘Š®†CÊnŸïÓœ÷Sáß" },
	{ "ya",L"ñâÞë›Å’¥Ý‘Q«ežõÑ¾åEùsÑÂý\’~—âˆR¬ˆ‚oˆº—¿ë²†s¶–ÑÅèâˆ×Ž†ˆÛªcÒˆBÑÄélÑ¼¸E‹I†¡¯{æ«Ûëþ†Ñ¿Š´Ñº…ìø†ù“Ñ¹øfÑÆÔþ…’ŒS„²çðÑÁâX…|á¬Ñ½‰ºÑÇªmšå¯PˆLçŒÑáÜˆÑÃý…–‘ùg„Ê‹’é…ƒÑ» ëØóŽÞíýÓ ÑÀÑÈÂyŽâåÂè›ðé·Ší¼¸Ž" },
	{ "yai",L"ñâÞë›Å’¥Ý‘Q«ežõÑ¾åEùsÑÂý\’~—âˆR¬ˆ‚oˆº—¿ë²†s¶–ÑÅèâˆ×Ž†ˆÛªcÒˆBÑÄélÑ¼¸E‹I†¡¯{æ«Ûëþ†Ñ¿Š´Ñº…ìø†ù“Ñ¹øfÑÆÔþ…’ŒS„²çðÑÁâX…|á¬Ñ½‰ºÑÇªmšå¯PˆLçŒÑáÜˆÑÃý…–‘ùg„Ê‹’é…ƒÑ» ëØóŽÞíýÓ ÑÀÑÈÂyŽâåÂè›ðé·Ší¼¸Ž" },
	{ "yan",L"œÑÙ÷Ðöoú`î»Ú¥ƒ°mŽtÑÖ¿t‡ÀãÆ‰ÁÑ×©éŽŸÌŠ°ÓƒÑÕáD…]ÙÈ™LáÃÑÓ…yî›ýzøeüsÑå†Íª_âûŒEéZŸSòžÒÚI«ŠŸeôeßVØWŽrµhŸðÑà³x‰üüiž¹äÎƒ¼òzÔPîƒüjýŒÑâØßô|äÙÓ_¼³šŠÔåûŽMý‡ÑÜÑæîÆGÖVé‘ÒóÈŠêšœÄÇ¦ÑßÑsÑÛÌš÷úÑËüfš‡ëÙ•à‰ÌÆFóÛÜ¾Žv—¦ÙðŸŸàIá€‘±‡²×…ªP¬JÑäþŸž Š×ìÍý]”©ýdÛ³ÑéÓ…ûšû}ž·ük‰†ÑÏéÜ½žmÑÔÑØØÍ‘î™¿°‹éÔÙž—ãázøHúŽÑç®[ÑÝ wæÌÑáÜyÑÍé’ZÑÞÙ²œ{üdÎiƒBØÉÝÎÑÚžÏìv“CÑÐ•VŠ¶î††«ÈTãÕ’ïž¥‚¹ ²ãUv³ŽÅE…——â×—ÑÉ†Çº™»Û±þÑÎëC›¡ÚÝýBŽi›WÁw›þÈCÑÒÑÊêÌ÷Ê…’…˜óFà‡{òVõ¦•óŠzðÑèëçŽs‘Ã‹Ç—ðÇráZùž‹jºcçüŒß“RÈ€öØVÊB‹÷ÑŠø‘û’´NÑÑ‘þ•¶ÚçÑãá‰ÑÌ‚©™•Ý" },
	{ "yang",L"ãóÑñê–Ôhï…Ñó§ÑîÁk±j«`Ñ÷ø—‘ÄöuÃo¶@¯ƒ˜Ó˜”è–¤Î^å}žYÁf«ŒÖ„½ÑòãZµSšç”aïrÑöÑìÑð‚ê’tÑë—î–³êgÑõˆ””®ñí¦ÑúâóûF¬„Œ÷‹Pžæ…n°W•DáàìRÝIë‡ÑïÑôšTûðB˜DÑùø„ï^ë›÷±“PšÞ…ó½D„ØòÕì¾ˆtÑíéAŠšÖUŸ¬ç{ìÈÑêåÑø•ª•[±ˆÝŒ" },
	{ "yao",L"÷¥è€Ò§š¥É@ËŽÌi¸GÑþ‚çø€ÑüœÈŒaæcŸÆÛuŒ¸ ú¹OÒ¦çÛÅ±â_Ö{ÙÒ¤¯‘‚¶úrœøðPª’Ò¡ÑûÈ™×Š²‡†ºÒ©Ñý•ê“ÁáÊ¬ŽÌÕµn–ÌŽAôíÒ«ƒeÔ¿ÚŒÊçðÎÀf•¬Ô@Ð‰ÝU dï¢Ø²Ôo‚xÒ¥“êŽCÆwëÈöŽ‹„’qé™Ëa˜l½Ä“eÒ¢ø^ïuê×ªqçòßb·šžìŒëýoèÃÒßºî–ž÷iÒ£ïŸòˆÜé‹Q†ÚÖ|ªrï_˜eˆòÒªØ³Ô¼É|«Q¼s¸Hì‰“uš|„üýGáæ·Žfò[·—ê´t–”ˆ±lÃ´ã“Ò¨¦ñºé÷Å—ú_" },
	{ "ye",L"‰¢”@ìÇæEí““ü {ügl•Ïùw£Ò¸Ò³ð†˜Iîô‡S‡™ÞÞ•¢š†­LŸ¤Ò±›˜G²|Ò¬•öˆìäyæU‰­È~ŽI’wñ@þ†èHÐJÚËŸºðvµB”IÑÊ†œŽJ ”ÚþâXÒ²ædÒºÖ]‚œ–¥Ò°²wà’‹­’šS–‘ƒp×§êÊÒ­–¦cÒ·š‡ûE•â‚´…½ŸîÐMÒ¶”KØÌóBÒ´Òµ°‡ðYÖÐ°”LìvÒ®àvÒ¹’Å’À¯u•ÐÒ¯ˆ¸" },
	{ "yi",L"îUÞÄÑ‹ÛüÌš…È^Þ²ãŽàæ‘üß½‹f‰©ÈU™}ð†‚XŸÁÅœÒÄ’íq@üpàÉÒâŸÛ…êöGøCÉšUÒÎ°¬ WŒT‘÷›n öçËßÞÒ»Í†î‰ÒÆÛDÆiõlã¨ÙŒ–Ø{‰ßÒáÑ`“Ìæ–ªÁpÒÞ–¤ß×ù€ºmÞ~Âk”¾ÙOÒÊŒh—×ÒÚÒ¿¿ˆ÷ðË‡ÜÓÆNÐ„Ø\ì½á{¥áÚÒÁÐ‘žË¶hÍ~¯ŽÎœúsôýátìˆ›¶õkŸ|™Ž•Æqñ´Ñvàc—©–õãžÒËôà›ÅÒßÒîÒëá»ñ¯ÒÔ[ý~Î’ÒæÔ„îÆÒäïîÁrÞTÚ˜ô¯Ê³åÆÞÚÜèêÝƒ|½X„ÖÒÍÔm•ËðêÒÕÒ×±Ž¯œ™ÒàâÂÒÅ J²e‡ÒìÚ’Ø[–ŽFÞÈØ×«p×g•âùÒÏìJ”î‚ëc–sÔrˆ`ˆË…À ô—àƒÞêe­Cƒx”ª~Ó”™ýM›¥çFŒ•âNŸÖîÐš™öŽå·Fµtš]…F…å›ª„·ŸyôèƒŒƒÏÛÝéìVÏï±Œb®A˜¯û@š¡•iŽƒÎõÞ æ„óAŽ–éóÒÉÜ²ß®–pÒÐ˜]Ö–Ù«ÒÓâøo£ÒA„ãÀ[–ñÀXèŸúœìaË„šcÒã‹ÎÐz•”ãAÒÀ›ÎÅ’Î•áyâzîVÝWø˜Š‰Òé‹ÂÒ~²G¤Ò¼‘«Äjï×ˆIûkÒÒÌ[ñÂü]ÕxÒÇÒÛ¯–Òèá¸”›uµKãi«}Ò½‘›ÒçâPûpðùÛÒÙØs×” D®×bš­í›ÝrÒïÒêØý†jßzûoÎ²©×‚‹¡ÔTÒÈÒìârù•Ù–åÐtÒÝòæØŠ˜à¹•Â]Ò]â¢ÒåÍ‚äô‹Òíú…Þ–ÃEÒ¾ÞÒÃåW™jÝ}ØæÚ±úgÒÖÓ~íô”§×hž‹·jÒÑú^Ø—î{¸v‹Úó`Ù“æä…¥ØîÁxÒÂêdÞjœjê‹ÒØµEÊÌˆ“~ŸéÕB¿O¯m‘ýÒÜÒÌèO”¹ýt–¶BÉß’L‚Ãˆ£×r¶êŽKŒ–" },
	{ "yin",L"¶†µš‡‘›ÆgÒòð°EœšÒõì‹AóS¹Nì‚áS|›Žà³þ“‘@ ì¿tõg™ƒˆøêŽ°aê”†‚ä¦½s‘ÓÛÈÒùãy™’î÷Ü§Ò÷âYÄŒíÒðâwñ«—VÖNÔCÓ¡Ó—‘\ÒöáDÒø‹HµåµœÞ´€ˆ¤âiö¸‡¨ªZÌa×ãŸ‡àšP³wÚ_ö¯ï‹ÏPÖ•Òú­KÏrýlë[ž@™ÓÒýÚyê›œôØ·ï‡ŸÝl–@¡äÎÑPÒûò¾]ñ—‘¶™aë³Òüâ¹úžôéž‘€ý‡ƒÜë–Òñ¯ŠŽ\Òó´Hë ÊaëLñ¿‡ôÜáÊ_ZþÛßÒôÕzœ^J¾žáþƒøÛ´Ûóêf™ýÉM‡wÒþÇZ”Õš’«lßÅ–ðÑÌ«ý]ˆŠ" },
	{ "ying",L"íŒ Iž„ÓLævÙa”lõ­‹‹ëúD¬“ƒOÓ¢™Õ–P‰à×GŸ–À†ñ¨Û«›sœ»À›Ó­Îží‹Žc›Æ†¦·f»Y_Ó£•£´Q‡Âœîì™ûWÓ¦Ó°îeÌcðÐ‚\éAïI•@×sž‰öäÞžLÜãžuÙøŸÉçø‚Ÿ®OéºúLž]‹””tž­Ètµ_ÓªÓ¨è]çÓ±°`è¬ÚA}ÝÓâßŒ[®Zå}“²ÖhŸ‡ò£äëÓ«ÐNûK÷j‡|ÞüÝº¿MÀtÜþÓ²Ó®—HÓ¬ëô‹kÓ¥‹ýúˆÝöàÓ‘ªŽgå­ú—°Ÿœ€Ó©Ñš—w»kÎs™ÑÓ¯ó¿ë›¾x‰LÓ¤œÁÓ§ÊËpÀ”Ï‰Ó³³AÄ{" },
	{ "yo",L"à¡‡©Ó´Óý†Ñ" },
	{ "yong",L"ÛÕÓ¸‡‡÷«à¯Ôò„ÊbÓ¼ð®“í¾‰M÷Ià{úçß†Þ°bÓÃÓÂÓÀÓ½ÛxÓ»žœÙ¸÷‘˜Ÿã¼‚æÓºçOÜ­ïÞ¶HÓ¶öÓ¹Ó·ákúxÆoÓÓµëtïJÓÓ¿õ—“N œ³‹~³l¹c÷ÓäVÓÁ‹£­ˆ¬œ¥‚òÉKÑKî„–Ô°MÓ¾àa" },
	{ "you",L"œ±‘ÉÉ‰–ëû~ÓÎÝµ›ÁÍY®hÞœòøÁm™¢«DÓÒÆh ¶HÃ‘”åÓÐÓ×ÓÅ¶xªqß[—`šüMƒž‘î›|ßˆ ¨Ã…ÓÄ›Y™ÔÞÌÞ”žXÀläPÓÊÍœâ™QÓÓÓÕØüƒÜÃUà]òößÏÓÌðà¯_òÊÇx ûÓÆµvàóáRÝ¬ŠmÝ¯˜Aôœ†NÂuñfÊ~‘nôíå¶¼n÷îÁgà›ŠµÑ„ÓÍÙ§ÕT˜©ÓÑÂir†eµ™ŽîöÏë»ŒMJŒ‚ºòÄ÷†õOÓÇ“AîðÔIféàäBJ—XÝjßKÝ’ØzèÖÁhÓÏØÕ‡¦ÓÔÓÉÓÖÓÈ÷øÓË" },
	{ "yu",L"–ëÓÞ‹žðNÓ÷µ€ÓëáCÆ‘îÚÓî…PÓèÓñßŽ„žÁ¹z—åƒÎµœùì¶·YâÅÁN°K@—§˜K²œÓÚ™óÉfO»BâDƒh«]–üô~ÖIÓð¾sôr÷Nï„É™ˆïØñ˜@êÅø\å÷òÊ ø…·{ÓÜÓãªz”ùòe†¸ý{­ÎC ŒókŸ~ÓûçŸ¬Z‚Rö§ä`Ü†Óú»nÈhÓßí²ßNìMìÐÔ¤ÅcäoŠÊÞíôcðõÓØþ–±ÍGâÀõ‚å[£ð|òâ¿›Ý›ÎƒöiáÎÕ˜é‘îAæúê|ƒÊ³‘þ’™ä´›šuý›Ô¢Î¾Ëv ¢ÝÒˆèšeãƒÓÝïJðÖå“Ø‹¬^¬r¶rùOóÄÄøˆû‡ÓàÖàÓÙÏLì£ûC†‰‘jžºìÙîY‚qàNî„œMÓæÊš…°Ž÷†ÉëéÅ„‘µšQÌ]‚øÝÇ—šØ…ŸúˆSÓç®Œ…ÇòõêìÓâÓþàôøƒÛu‹äÌPÊv–fºh‹VØ®ª—™Ñˆ³†³_·UÔ£ŽZ’GëkÆœÓüÓìèžðÁËØ¹ÆRó^ìÏ±EÓéµN®ÈgÓù¯ƒ™ÕZâ•¹ÈÓòÝ÷îÓåôdìÛú–åýâ×Ô¥‘íßyÁ|é“Óö}ÓÛÅ†›@Ù¶ÓóýrÞX†³·Ýh VÔ¦ö¹ç~·Cô§à¯ó‡oæ¥Óï·‹áqàhëTÓä‡‰ðöÑ@áüÚÍ”ËÓêÓDÁ”ñSûOÞ}­m»Z¼uÓø” ô¨ÚÄ”Ñ‹àö×uŠØêœãÐ‰¥÷r‚¦ÓíÓõÖ~ñÁˆÖ²IÓá¶RÐsÓô›Aú}€’§œŸœUè¤«_Œ†ÏXôˆÓý™Èñ¾’TÞzÔ¡" },
	{ "yuan",L"ÉVÔ°ò{ä‘Ø’ÎmîŠ‹õœ®Š††TÔ­Ô¹ÎzÐcŠ€ß‡Ô²±\Ôªæ…›ðß–Ô©Ô¯ƒÒßR…øSÔ¶ŒwÖwÑj€œmà÷­Þ@ˆA¸–zÔ«Ûùë¼ÔµÔ§ð°‡…‹…ÆŠ˜rËeÃO‡ûÔ·úM—¥Í›ªxûg™´øx„uœeù ö½Á~Ô¨è¥óîüŒ…ŒÔ¬Ü«˜gÔ±ÚOœaô’Ô®ËQùt‡äÔºªjÔ³éÚˆ@ž”ó¢‹…™œÆÎQ†¿¿Fü‚ÓÉdæÂÍWÜ¾ßhÑ†Ô¸˜Cíóüxñr‹ œYµž¾‰ãäÔ´ÉA‰íÑ“ÞòÑr" },
	{ "yue",L"Šxè€‚ügé‡ÍR–†Õh¹–ó–¶^ÙßßÜ™µîáÕf‘à»l•õšõãX ~«hÔÂËµÔÃ³E»Câ_ÌgûN§éÐé†ÔÄûVºMÍQÔ¼µjÔÀ˜·‹í¼s»aÔ¿Ü‹ÚŒÔÁ²ˆŒéÔ½Ô¾Ž[ÜSå®˜SÔ»‡‚ÄŸÚ”èÝÜVë¾’Õ’`x»›†dÀÖ¦" },
	{ "yun",L"ÔÆÊ•Ø’›éâqÔÍájÔÈë…ÀˆŠ@†TíMËœàyáñÔËýq˜øÉl¹ošèÉCÔÎëEÙšã³Ç\›VÌNÝ˜ìB¼‹˜XŠ[ä]ÑŽÎ‚•žírða†½»±dí‡çÛ©ñNóÞÙ„ÔÊÄZàiÔÌÔÇŸÂÝœêÀŸ±ÞdÉQ‚ÖœÝ‘Cýy«jÔÐä·¶nç¡ß\ëµÚSÊ|êméæã¢ŠuÔ±Ÿ¾ÚOÔÉÁÊŸè¹®sÜ¿ád¿Z‹‹ÔÅÂmíy„ò›âÎQÔÏîf–—ñaŸ¸’d¿aìÙ¾´pºJšŒ¿A" },
	{ "za",L"¼’›eëj›jÅNÔúÔÛ‚ÌÞ†¹ÅHŽ‰–ýßÆ‡mÔÒô˜‡Ôë{‡ÙësãNÒS´’ÔÑíˆ¼™Õ¦ÞÙÔÓ" },
	{ "zai",L"áPÔÚÔÕ¿fœ…ÔÙžüœÖ’Dƒ„Ô×çÞÝdÇáÌÔÔÔØÔÖ²P›’×ÐžÄÙ†‚î" },
	{ "zan",L"à™¶`çZžUÛŠÙ“SÔÞ“Ëè¶ÔÛºð•ƒ›‚Ì†¹­‘ƒ³ºdÙmç‘–ý­ç‡áA·‰‡ÔŒv×“”€çYÒ{öÉêÃô¢ÔÝàŸôØÔÜ×{ÚŽÞÙ•ºôõž£‡kƒ­" },
	{ "zang",L"ÔàÁn…MÙjn™âÞÊê°óG²ØÔßÙ_óvÚNñzÊiÔá‰Z ™ÅKÚEäQæà" },
	{ "zao",L"ÔâÔéÖ¸^ËkásÛ›—_ßð——Ôî°oÔì×Y†×‚óÅ–Òè­bÀRÔí­FÔèŸ¯‘VÔêçØºrÔãÔçÔåÚ‹ÔæÔä†rÔëÔï" },
	{ "ze",L"ØŸ¶štÕ‹Ž¾’¾ØÓÖ‰†¨È[ºjÉê¾›gÔóÔò‰÷óÐÔðÒ]šò´Ÿ’kÔñïŽý`˜Áàý‹¨ŽÙ¡²àßõœÚÙ‘‚È‡K“ñóåÌEåÅ²žô·²cýv†‡ØÆ›z•WÕ¦Ï„tÂd°ƒœõûB" },
	{ "zei",L"Ôôöf÷Œ÷eÏŒÙ\‘å" },
	{ "zen",L"ÔõÚÚ×P×U‡×" },
	{ "zeng",L"êµ¿•´Œ™I×ÛçÕ¿fÙ›Ô÷Ôö‰ˆÔøï­ä{­QÖŸôiŸå•û¾C³DîÀÔùà‹" },
	{ "zha",L"¼’òÆßîƒÔŒo®hßå’€÷‡íCŸ¤žÁ“«“’Ûz×A‡ÍœÑ–ÅíÄél†Æà© £ám¹€ëÔý¼™“ƒé«Õ¤ÔüÕ©ýOõ~Ôþ÷þÔpÕ§‚¼Û‚Õ¨÷Ôûðl¹†ÔúâÇ×õ’Ÿ²éÆzÓuÜˆÍlÕ¥„ž°•åŽ–¼À¯ŠL°šýv×QÕ¡Õ¢¤˜ÏðäÕ£Þêß¸õW„‘…~Õ¦Âd’s" },
	{ "zhai",L"íÎØŸ˜z…‰ã²àÕ¬’Æ‚ÈµÔ“ñ”È’nÕ®¼ÀãSÕ­Ôð‚ùýS…~™yÕ¯ódÔñÕªñ©Õ«" },
	{ "zhan",L"ûrä‘ðë•Þø‹¶¬”ö—£®¾`–î’€”Ø´DÌ›°œÕ·²üò–˜öÌœÍtì¹Õ²Û…–énÕ¸±KïQÕ°ÔaŽEÚj÷gÇ•Õ½šÖ˜^Õ¹Õ´šØÓOÕ¼×Õ³×`Õ¿á\GÒf‡~™ÙÕ»ürÕºËUðeÝšÝu¬Wø@×dô}ÕÀÞJûDÕ±ßücÛ@ïãÚÞ‚·ò ðŒîÕ¶øÕ¾Õµïs" },
	{ "zhang",L"•ÀÕÂ÷Já¤ÕÅÕÃØëÕÈÉŸÕÇð\Ÿâ¯´˜ÕÏqáÖÕÍ³¤» éL›îÕÊ{Ù~æÑ’Eè°Ž¤ÕÌû–¯“ÕÁˆ²dÕËÕÄŽÇÛµò†ÕÆßl‰zÕÎó¯‘Pƒ@Ã›¯oéMì ÕÉ" },
	{ "zhao",L"¸SÕÒŽ‚á“îÈÇŸ†Ô™˜ÖšÁ^³¯°œ  YÕÖ±@ªD×ÅÚwñq²ÕÓõeÚ¯žÝÔtèþÃAŠ„×¦ˆå™¹|¬ÕÑÕØ”í³°––ÕÕßúÕÐóÉøJÕ×ÕÔ•×ÕÙÃDü{ü…ãD" },
	{ "zhe",L"èÏñÞžâºÖøÝtúpéüÕáÕÝäO‡¬Ým×Å’VÕß^ÍE•‡ÂzíÝ×yô÷‹«õ„ÏVÕÚ–l˜Î»q•†æN³YÕÛ×„†Ö†ÐŸ†øÞHÖ•òØÏUðÑ³K‘eúvÒxóCµ—Õã…zšy†£—‘ÕÜ†´ß@Õâ»„Ô€œJÚØß¡Õàˆ³ÕÞ" },
	{ "zhei",L"ß@Õâ" },
	{ "zhen",L"áGŒÇÕð•Håg±‡Èœ÷yÕó˜ˆ‚E’ræä‹–bæ‚–žÒ˜ÙcÉR±péô˜^Ñ]«‚Õç“Žé»•_ÀƒŽ¬Õò˜çÍ–’™´UÝFð¡‹‚É¼…ËmëÞÕèôIÕäÂrül“LœìkÕïåŒÕëóðá˜ÕåØé©ä¥âœøcêâÖnØ‘±wš‹ê˜EÕêÞtŠª”´çÇáIð²ÕæÛÚümµ½GÖ¡–×ü‡½„›ßZÚfÕgèåÐÕìÕéÝèÞÕñª€î³ŒzäÚ¬‘¿bëÓ‰`Ô\Õí›l¿j»Eê‡¸t—FæP¶GÝŸÕî–ÚŽžìõñ}" },
	{ "zheng",L"ï£’cŠ’ã`õ›Õô•“Õö˜ÍÕûÕú’êàóÝñÕøÔ^’ð”Ö¤šéõSŸAåPÕùˆÁÛt‹oîÛ“Õ¶¡¹~Ö£±kˆ½Ñ Ž‰^áçÕþÚºÂtºPþÖ¢± ‘~ÕŠøg»¼lþÕýòÕ÷ÃwÕüžÚöëÕá¿ô@°YÕõªb“@×C" },
	{ "zhi",L"’†ÃŽæÖÀ…æšl˜Þ¶~¹eÜW‘ç’ÃÂpËŒÖ¸Ð}Ž‚À†¶ƒŒ…ïô ÃÖ­…„›b®‡¶q‚Ž¶ˆõÅD—Ð”òòÖ®áçðëèv‘ÁªOâåÛ¤½’WÇ õôÆW—„ÖÄÖ°Ö¬ªöSÖ½ÖÁ¶héòŠ©Ïdé@„¶¼•Ö»ÜUµ•ÛNÃqŽŽ±‚Æ‡ÄˆÖÌÙ—–»Ö¶âºÓd›±Ìuœ]ÖÎÐ—Û—íéÛ•èäÖ¦Ø´Ö¯ã‡‘Æ½‚µwÖÃ”`Ö±Ö¼ÖÍæïÜÆŽèÞý·W¶_ìíèÙ¯WÙ|Ñuìó„¬Ö«µ…–sÖÆÃeèeêeÖ¥åéúEúvÞŒÂšÖÈÖ¹Á“ÖËˆÌ”ÕÕñ\ÖÇžøT‘pÝeŒˆ€˜uòsÖÅ‚ÐÕI—ù¿@Öº‰yÏH„Œ–ñ™±èœ›œ¯FÍV¼ƒœZè×•y”SõÙ¶žž\‰~›E›‚ŸÜˆpÚìÛú”TÖ·¶oŽÃŠ‰’XÓhËSÔJ’nÖÂ¶AòÎÖ´ö£·aØTñcÖÊèÎœþÊÏÓzœFÐ˜Ö©ÜÖ³õÜ—d“¯ŠÍåëdàùëÕ¼ˆü~ÖÏ™£ëb­•Å]õ¥ðºÆÖ¿ëù“w¿{ÝTµYÖµ˜àòŽØ ßtñ‹ên’”Ö§‚uØUÖ¾ÖÉ“ˆÊ¶ÒžÒjá™³UÖ¨¯€œí×RäkÖ²êÞ‹Àˆ^¯U‘eôêäK¿—øvéù­}ÖªªaÅ\“´ø›DµoÐ" },
	{ "zhong",L"Šq«ÖÒ³ïñ¾…ŸŽçŠÖÕÊW·NÐxÐ\ÏxšpÖÙÖÖ‰V·rÖÓÎuÖÔ–°Í\ÖA ðÚ£»bÖÐ›Oµr†ÁŠtâ{›wÄ[ˆúÛ Ö×æRó®‹gô±×¯~½KÖÑ±ŠÆ žÆŽºäVâìø‚õàÖØ„d‚£Î ÖÚŒ»Ð{¹Wü™â`" },
	{ "zhou",L"±TÖæÖèòLÝS°™¿UÖå…âþVóEÖÛôüæ¨Ý§æûÝcÁŸ×žÆ®L‡œëÐþ_×pÕŒ‡€—¹œ@ã{þUù@ÖÞÚQ¹ÈF»Q³B»‹Ëg›öBÝqžëÖáÖâÅ«‰íØ†Bú•ŽƒÙôíûbÖÜÙk•ƒžöë“ç§þ`¼qÖçÈ’ñt†µÖàÖäô¦ÖßßúÖãàXßL‹B»N²HÞb¯JÖaÖÝÔkÔ—ƒu" },
	{ "zhu",L"ô¶„±•ôÊôÀ‚™½äóÇAÇdôãÁqñ[Ùª”±ÖéõfÎwË\–XÖúŒFÖTðñ¼Ÿä¾äŠÐEñvÛ¥×¡èTÚŸèÌ¹hÆ^Öì”²ê• ‰žÛ²šè“ò|ü}éÍóç¸˜ÝO¶‹–Çðæ·”ŒÙµ‚Öî„¸³p×£±vÖò½AÐÛHÖñ÷E×¢›{™Á½ZúžÊõÖøŒeñÒãLÙAÖö³d”áÖýÏŽÜïË ÷æñ–ãËŸÖêÜÑ™·ˆ|ÆrÐWë—ÑNÖ÷šŸÖëßI–’ö^ÕDÖþÊxéÆÖíä¨î¸mŒ¥Ô}ØùÖôÉ T‰Ô­ÁCºaîùÔ]¸‰×¤ÖüóÃž¯ïŒø–Ÿ—‡ÚžzÖûÖó™îìÄ˜ÖÞŽÖðõîû„ÖïÖõºZØiÖùºB" },
	{ "zhua",L"×¥Ä“ëÎÎºœ™tó˜“«×¦" },
	{ "zhuai",L"ÛJ×§±‘àÜÞD’Åî“×ª" },
	{ "zhuan",L"ÂZ`ÞDƒQÖKÉEžÀ»Mãç×«ò§×¬þz„–‰tŒNâÍ®UÙ´s‡Ê×­ƒ]´u×ª¬ƒ…¡÷HÒNº‹Ïm¿x×©×¨àî…ºeÄxºiˆæ´«‚÷ð‚ÄR­A¸|‹§Œ£×Nßù" },
	{ "zhuang",L"ÞÊŠÏÑb×±×°×®ŽáŠy—[»’˜¶‰öœ³Çf×¯¼P" }
};