/*
 * Offline participant component-report viewer.
 *
 * Usage: build_ocular_component_reports.py embeds this source and REPORT_DATA
 * into one participant HTML. Open that HTML directly in a modern browser.
 * No server, network request, Python session or scientific file is needed.
 *
 * Data flow: select a run -> decode its compact float32 peak samples/uint16
 * native offsets on demand -> draw only the visible time window on canvas.
 * Times remain seconds from the exact recording start. NaN bins are explicit
 * unscored gaps and never connected. The displayed values are normalized,
 * peak-preserving diagnostic copies; detection always uses original caches.
 * Browser controls change display only and never save scientific decisions.
 */
"use strict";
(function () {
  const data = JSON.parse(document.getElementById("report-data").textContent);
  const el = function (id) { return document.getElementById(id); };
  const ready = data.runs.filter(function (r) { return r.status === "ready"; });
  const decoded = new Map();
  let run = null, samples = null, start = 0, width = 20, gain = 1, drag = null;
  const canvas = el("traces"), ctx = canvas.getContext("2d");

  // Decode only visited runs; Float32/Uint16 values are little endian.
  function binary(text) {
    const bytes = atob(text), out = new Uint8Array(bytes.length);
    for (let i = 0; i < bytes.length; i += 1) out[i] = bytes.charCodeAt(i);
    return out.buffer;
  }
  function decode(r) {
    if (!decoded.has(r.source_filename)) {
      const values = new Float32Array(binary(r.display.values));
      const offsets = new Uint16Array(binary(r.display.offsets));
      if (values.length !== r.n_components * r.display.n_bins || offsets.length !== values.length)
        throw new Error("Trace payload dimensions disagree with the report.");
      decoded.set(r.source_filename, {values: values, offsets: offsets});
    }
    return decoded.get(r.source_filename);
  }
  function clamp() {
    width = Math.max(Math.min(.25, run.duration_sec), Math.min(width, run.duration_sec));
    start = Math.max(0, Math.min(start, Math.max(0, run.duration_sec - width)));
  }
  function controls() {
    clamp();
    el("time-slider").max = Math.max(0, run.duration_sec - width);
    el("time-slider").value = start;
    el("jump").value = start.toFixed(2);
    el("time-label").textContent = start.toFixed(2) + " to " +
      Math.min(run.duration_sec, start + width).toFixed(2) + " s";
    el("previous-run").disabled = ready.indexOf(run) === 0;
    el("next-run").disabled = ready.indexOf(run) === ready.length - 1;
  }

  // The canvas is redrawn from compact arrays, never from external files.
  function draw() {
    if (!run) return;
    controls();
    const ratio = Math.min(window.devicePixelRatio || 1, 2);
    const w = Math.max(360, canvas.parentElement.clientWidth), h = 36 + 64 * run.n_components;
    canvas.width = Math.round(w * ratio); canvas.height = Math.round(h * ratio);
    canvas.style.width = w + "px"; canvas.style.height = h + "px";
    ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    ctx.fillStyle = "#ffffff"; ctx.fillRect(0, 0, w, h);
    const left = 100, right = w - 18, scale = (right - left) / width;
    ctx.font = "12px system-ui"; ctx.textBaseline = "middle";
    for (let tick = 0; tick <= 5; tick += 1) {
      const x = left + tick * (right - left) / 5;
      ctx.strokeStyle = "#e7ebf0"; ctx.beginPath(); ctx.moveTo(x, 30); ctx.lineTo(x, h); ctx.stroke();
      ctx.fillStyle = "#536174"; ctx.textAlign = tick === 5 ? "right" : "left";
      ctx.fillText((start + tick * width / 5).toFixed(width < 2 ? 2 : 1) + " s", x, 15);
    }
    const step = run.display.step, sf = run.sfreq, count = run.display.n_bins;
    const first = Math.max(0, Math.floor(start * sf / step) - 1);
    const last = Math.min(count, Math.ceil((start + width) * sf / step) + 1);
    for (let c = 0; c < run.n_components; c += 1) {
      const top = 33 + 64 * c, middle = top + 31, base = c * count;
      ctx.textAlign = "left"; ctx.fillStyle = "#172235";
      ctx.fillText("ICA" + String(c).padStart(3, "0") + "  [" + c + "]", 8, middle);
      ctx.strokeStyle = "#edf0f4"; ctx.beginPath(); ctx.moveTo(left, middle); ctx.lineTo(right, middle); ctx.stroke();
      ctx.save(); ctx.beginPath(); ctx.rect(left, top, right - left, 62); ctx.clip();
      // NaNs represent entire unscored bins, including gaps shorter than a bin.
      let gap = null, line = false;
      ctx.strokeStyle = "#172235"; ctx.lineWidth = .9; ctx.beginPath();
      for (let i = first; i < last; i += 1) {
        const value = samples.values[base + i];
        const time = (i * step + samples.offsets[base + i]) / sf;
        if (!Number.isFinite(value)) { line = false; if (gap === null) gap = i; continue; }
        if (gap !== null) {
          ctx.fillStyle = "#e6e9ed";
          ctx.fillRect(left + (gap * step / sf - start) * scale, top,
            (i - gap) * step / sf * scale, 62);
          gap = null;
        }
        const x = left + (time - start) * scale, y = middle - value * 26 * gain;
        if (line) ctx.lineTo(x, y); else ctx.moveTo(x, y);
        line = true;
      }
      if (gap !== null) {
        ctx.fillStyle = "#e6e9ed";
        ctx.fillRect(left + (gap * step / sf - start) * scale, top, (last - gap) * step / sf * scale, 62);
      }
      ctx.stroke(); ctx.restore();
      ctx.strokeStyle = "#d9dfe8"; ctx.beginPath(); ctx.moveTo(0, top + 63); ctx.lineTo(w, top + 63); ctx.stroke();
    }
  }
  function setMap() {
    const kind = el("sensor-kind").value;
    el("maps").src = run.maps[kind];
    el("map-caption").textContent = kind === "mag" ? "Magnetometer component maps" : "Gradiometer component maps";
  }
  function selectRun(source) {
    run = ready.find(function (r) { return r.source_filename === source; });
    if (!run) return;
    samples = decode(run); start = 0; width = Math.min(20, run.duration_sec);
    el("window").value = run.duration_sec < 20 ? "full" : "20";
    el("run-select").value = source;
    el("source").textContent = run.source_filename;
    el("run-summary").textContent = run.n_components + " components · " +
      run.duration_sec.toFixed(1) + " s · " + run.display.display_rate_hz.toFixed(0) +
      " display bins/s · convergence: " + String(run.converged);
    el("sensor-kind").replaceChildren();
    Object.keys(run.maps).forEach(function (kind) {
      const opt = document.createElement("option"); opt.value = kind;
      opt.textContent = kind === "mag" ? "Magnetometers" : "Gradiometers";
      el("sensor-kind").appendChild(opt);
    });
    el("note-row").value = data.participant + " R" + String(run.acquisition_run).padStart(2, "0") +
      " | blink: [indices or NONE] | eye: [indices or NONE] | notes: ";
    setMap(); draw();
  }
  data.runs.forEach(function (r) {
    const opt = document.createElement("option");
    opt.value = r.source_filename;
    opt.textContent = "Run " + r.acquisition_run + (r.status === "ready" ? "" : " — " + r.status.replaceAll("_", " "));
    opt.disabled = r.status !== "ready"; el("run-select").appendChild(opt);
    const tr = document.createElement("tr");
    [String(r.acquisition_run), r.source_filename, r.status.replaceAll("_", " "), r.summary || ""].forEach(function (text) {
      const td = document.createElement("td"); td.textContent = text; tr.appendChild(td);
    }); el("coverage").appendChild(tr);
  });
  el("participant").textContent = data.participant + " · Ocular components";
  el("report-id").textContent = data.report_id;
  el("run-select").onchange = function () { selectRun(this.value); };
  el("sensor-kind").onchange = setMap;
  el("previous-run").onclick = function () { selectRun(ready[ready.indexOf(run) - 1].source_filename); };
  el("next-run").onclick = function () { selectRun(ready[ready.indexOf(run) + 1].source_filename); };
  el("time-slider").oninput = function () { start = Number(this.value); draw(); };
  el("jump").onchange = function () { start = Number(this.value) || 0; draw(); };
  el("window").onchange = function () { width = this.value === "full" ? run.duration_sec : Number(this.value); draw(); };
  el("gain").onchange = function () { gain = Number(this.value); draw(); };
  el("back").onclick = function () { start -= width * .8; draw(); };
  el("forward").onclick = function () { start += width * .8; draw(); };
  function zoom(factor) { const center = start + width / 2; width *= factor; start = center - width / 2; draw(); }
  el("zoom-in").onclick = function () { zoom(.5); };
  el("zoom-out").onclick = function () { zoom(2); };
  canvas.onpointerdown = function (e) {
    drag = {x: e.clientX, start: start}; canvas.setPointerCapture(e.pointerId);
  };
  canvas.onpointermove = function (e) {
    if (drag) { start = drag.start - (e.clientX - drag.x) * width / Math.max(1, canvas.clientWidth - 118); draw(); }
  };
  canvas.onpointerup = canvas.onpointercancel = function () { drag = null; };
  canvas.addEventListener("wheel", function (e) {
    if (e.ctrlKey || e.shiftKey) { e.preventDefault(); zoom(e.deltaY > 0 ? 1.25 : .8); }
  }, {passive: false});
  window.addEventListener("resize", draw);
  window.addEventListener("keydown", function (e) {
    if (/INPUT|TEXTAREA|SELECT/.test(document.activeElement.tagName) || !run) return;
    if (e.key === "ArrowRight") { e.preventDefault(); start += width * .8; draw(); }
    if (e.key === "ArrowLeft") { e.preventDefault(); start -= width * .8; draw(); }
  });
  if (ready.length) selectRun(ready[0].source_filename);
  else { el("review").hidden = true; el("empty").hidden = false; }
  document.documentElement.dataset.reportReady = "true";
})();
